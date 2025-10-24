/*********************************************************************************
 * MIT License *
 *                                                                               *
 * Copyright (c) 2025 Sang Chul Choi & ChatGPT * Based on Oatk code structure
 *(c) 2022 Chenxi Zhou                             *
 *********************************************************************************/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#include "ketopt.h" /* CLI parsing */
#include "kseq.h"   /* FASTA/FASTQ reader for gzFile */
KSEQ_INIT(gzFile, gzread)

#ifndef SYNCFILTER_VERSION
#define SYNCFILTER_VERSION "0.4.0"
#endif

/* ============================== Small Utils =================================
 */

static void die(const char *msg) {
  fprintf(stderr, "[fatal] %s\n", msg);
  exit(EXIT_FAILURE);
}
static void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p)
    die("OOM");
  return p;
}
static char *xstrdup(const char *s) {
  size_t n = strlen(s) + 1;
  char *p = (char *)xmalloc(n);
  memcpy(p, s, n);
  return p;
}

static int makedir_p(const char *path) {
  /* create parent directories if needed (best-effort, POSIX-ish) */
  if (!path || !*path)
    return 0;
  char *dup = xstrdup(path);
  for (char *p = dup + 1; *p; ++p) {
    if (*p == '/') {
      *p = 0;
      mkdir(dup, 0755);
      *p = '/';
    }
  }
  mkdir(dup, 0755);
  free(dup);
  return 0;
}

/* nt -> 2-bit; -1 for non-ACGT */
static inline int nt4(char c) {
  switch (c) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  default:
    return -1;
  }
}
static inline uint64_t rc2b(uint64_t x, int l) {
  x = ~x;
  uint64_t y = 0;
  for (int i = 0; i < l; i++) {
    y = (y << 2) | (x & 3u);
    x >>= 2;
  }
  return y;
}

/* ============================ Closed Syncmers ================================
 */

typedef struct {
  int k, s;
} sm_param_t;

typedef struct {
  uint64_t *a;
  size_t n, m;
} u64v_t;
static void u64v_init(u64v_t *v) {
  v->a = NULL;
  v->n = v->m = 0;
}
static void u64v_free(u64v_t *v) {
  free(v->a);
  v->a = NULL;
  v->n = v->m = 0;
}
static void u64v_push(u64v_t *v, uint64_t x) {
  if (v->n == v->m) {
    v->m = v->m ? v->m << 1 : 1024;
    v->a = (uint64_t *)realloc(v->a, v->m * sizeof(uint64_t));
    if (!v->a)
      die("OOM");
  }
  v->a[v->n++] = x;
}
static int cmp_u64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}
static int cmp_u32(const void *a, const void *b) {
  uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
  return (x > y) - (x < y);
}

static void collect_closed_syncmer_s(const char *seq, int len,
                                     const sm_param_t *P, u64v_t *out) {
  const int k = P->k, s = P->s;
  if (len < k)
    return;
  const int ns = len - s + 1, w = k - s + 1;
  uint64_t f = 0, r = 0, mask = (s == 32 ? ~0ULL : ((1ULL << (2 * s)) - 1ULL));
  int l = 0;
  uint64_t *sh = (uint64_t *)xmalloc((size_t)ns * sizeof(uint64_t));
  int nw = 0;
  for (int i = 0; i < len; i++) {
    int b = nt4(seq[i]);
    if (b < 0) {
      l = 0;
      f = r = 0;
      continue;
    }
    f = ((f << 2) | (uint64_t)b) & mask;
    r = (rc2b((uint64_t)b, 1) << (2 * (s - 1))) | (r >> 2);
    if (l < s) {
      l++;
      if (l < s)
        continue;
    }
    uint64_t can = (f < r) ? f : r;
    if (nw < ns)
      sh[nw++] = can;
  }
  if (nw > ns)
    nw = ns;
  for (int i = 0; i <= len - k; i++) {
    uint64_t min = UINT64_MAX;
    int arg = -1;
    for (int j = 0; j < w; j++) {
      uint64_t v = sh[i + j];
      if (v < min) {
        min = v;
        arg = j;
      }
    }
    if (arg == 0 || arg == w - 1)
      u64v_push(out, min);
  }
  free(sh);
}

/* ============================= Folds / Classes ==============================
 */

typedef struct {
  double f1, f2;
  int three;
} folds_t;
/* Return: if three bins: 0=nuclear, 1=mito, 2=plastid; else
 * 0=nuclear,2=organelle */
static inline int classify(double x, const folds_t *F) {
  if (F->three) {
    if (x < F->f1)
      return 0;
    if (x < F->f2)
      return 1;
    return 2;
  }
  return (x < F->f1) ? 0 : 2;
}
static inline const char *class_str(int cls, int three) {
  if (three)
    return cls == 0 ? "nuclear" : (cls == 1 ? "mito" : "plastid");
  return cls == 0 ? "nuclear" : "organelle";
}
static inline double dmean_u32(const uint32_t *v, size_t n) {
  if (!n)
    return 0.0;
  double s = 0;
  for (size_t i = 0; i < n; i++)
    s += (double)v[i];
  return s / (double)n;
}

/* ================================ GFA Parser ================================
 */

typedef struct {
  char *name;
  char *seq;
  double cov;
} gseg_t;
static int is_all_digits(const char *s) {
  if (!*s)
    return 0;
  for (; *s; ++s)
    if (*s < '0' || *s > '9')
      return 0;
  return 1;
}

static int parse_gfa_S_line(char *line, gseg_t *out) {
  int n = 0;
  char *tok[32];
  for (char *p = line; *p && n < 32;) {
    while (*p == '\t' || *p == ' ')
      ++p;
    if (!*p)
      break;
    tok[n++] = p;
    while (*p && *p != '\t' && *p != ' ')
      ++p;
    if (*p) {
      *p = 0;
      ++p;
    }
  }
  if (n < 3 || strcmp(tok[0], "S") != 0)
    return 0;
  int seqi = 2;
  if (is_all_digits(tok[2]) && n >= 4)
    seqi = 3;
  out->name = xstrdup(tok[1]);
  out->seq = xstrdup(tok[seqi]);
  out->cov = NAN;
  for (int i = seqi + 1; i < n; i++) {
    if (strncmp(tok[i], "SC:i:", 5) == 0)
      out->cov = atof(tok[i] + 5);
    else if (strncmp(tok[i], "KC:i:", 5) == 0) {
      double kc = atof(tok[i] + 5), L = (double)strlen(out->seq);
      if (L > 0)
        out->cov = kc / L;
    }
  }
  if (!(out->cov == out->cov))
    out->cov = 1.0;
  return 1;
}

/* ======================= Lightweight maps (no khash) ========================
 */
/* u64->double coverage (graph modes) */
typedef struct {
  uint64_t key;
  double val;
} kv64d_t;
typedef struct {
  kv64d_t *a;
  size_t n, m;
} kvv64d_t;
static void kvv64d_init(kvv64d_t *v) {
  v->a = NULL;
  v->n = v->m = 0;
}
static void kvv64d_push(kvv64d_t *v, uint64_t k, double val) {
  if (v->n == v->m) {
    v->m = v->m ? v->m << 1 : 1024;
    v->a = (kv64d_t *)realloc(v->a, v->m * sizeof(kv64d_t));
    if (!v->a)
      die("OOM");
  }
  v->a[v->n].key = k;
  v->a[v->n].val = val;
  v->n++;
}
static int cmp_kv64d_key(const void *A, const void *B) {
  uint64_t x = ((const kv64d_t *)A)->key, y = ((const kv64d_t *)B)->key;
  return (x > y) - (x < y);
}
static void kvv64d_sort_unique_max(kvv64d_t *v) {
  if (!v->n)
    return;
  qsort(v->a, v->n, sizeof(kv64d_t), cmp_kv64d_key);
  size_t w = 0;
  for (size_t i = 0; i < v->n;) {
    uint64_t k = v->a[i].key;
    double best = v->a[i].val;
    size_t j = i + 1;
    while (j < v->n && v->a[j].key == k) {
      if (v->a[j].val > best)
        best = v->a[j].val;
      j++;
    }
    v->a[w].key = k;
    v->a[w].val = best;
    w++;
    i = j;
  }
  v->n = w;
}
static int kvv64d_get(const kvv64d_t *v, uint64_t key, double *out) {
  size_t L = 0, R = v->n;
  while (L < R) {
    size_t M = L + ((R - L) >> 1);
    uint64_t k = v->a[M].key;
    if (key == k) {
      *out = v->a[M].val;
      return 1;
    } else if (key < k)
      R = M;
    else
      L = M + 1;
  }
  return 0;
}
static kvv64d_t build_sm_cov_from_gfa(const char *gfa_fn, const sm_param_t *P) {
  (void)
      P; /* P is only used to compute syncmers; we recompute per S-seq below */
  FILE *fp = fopen(gfa_fn, "r");
  if (!fp) {
    perror(gfa_fn);
    die("cannot open GFA");
  }
  kvv64d_t M;
  kvv64d_init(&M);
  char *line = NULL;
  size_t cap = 0;
  while (getline(&line, &cap, fp) != -1) {
    if (line[0] != 'S')
      continue;
    gseg_t s;
    if (!parse_gfa_S_line(line, &s))
      continue;
    u64v_t v;
    u64v_init(&v);
    /* We *do not* HPC at GFA segment level — keep coverage as in graph */
    collect_closed_syncmer_s(s.seq, (int)strlen(s.seq), P, &v);
    if (v.n) {
      qsort(v.a, v.n, sizeof(uint64_t), cmp_u64);
      size_t m = 0;
      for (size_t i = 0; i < v.n; i++)
        if (i == 0 || v.a[i] != v.a[i - 1])
          v.a[m++] = v.a[i];
      for (size_t i = 0; i < m; i++)
        kvv64d_push(&M, v.a[i], s.cov);
    }
    u64v_free(&v);
    free(s.name);
    free(s.seq);
  }
  free(line);
  fclose(fp);
  kvv64d_sort_unique_max(&M);
  return M;
}

/* u64->u32 multiplicity (quickview) */
typedef struct {
  uint64_t key;
  uint32_t val;
} kv64u32_t;
typedef struct {
  kv64u32_t *a;
  size_t n, m;
} kvv64u32_t;
static void kvv64u32_init(kvv64u32_t *v) {
  v->a = NULL;
  v->n = v->m = 0;
}
static void kvv64u32_push(kvv64u32_t *v, uint64_t k, uint32_t val) {
  if (v->n == v->m) {
    v->m = v->m ? v->m << 1 : 2048;
    v->a = (kv64u32_t *)realloc(v->a, v->m * sizeof(kv64u32_t));
    if (!v->a)
      die("OOM");
  }
  v->a[v->n].key = k;
  v->a[v->n].val = val;
  v->n++;
}
static int cmp_kv64u32_key(const void *A, const void *B) {
  uint64_t x = ((const kv64u32_t *)A)->key, y = ((const kv64u32_t *)B)->key;
  return (x > y) - (x < y);
}
static void kvv64u32_sort_merge_counts(kvv64u32_t *v) {
  if (!v->n)
    return;
  qsort(v->a, v->n, sizeof(kv64u32_t), cmp_kv64u32_key);
  size_t w = 0;
  for (size_t i = 0; i < v->n;) {
    uint64_t k = v->a[i].key;
    uint64_t sum = v->a[i].val;
    size_t j = i + 1;
    while (j < v->n && v->a[j].key == k) {
      sum += v->a[j].val;
      j++;
    }
    v->a[w].key = k;
    v->a[w].val = (sum > UINT32_MAX ? UINT32_MAX : (uint32_t)sum);
    w++;
    i = j;
  }
  v->n = w;
}
static uint32_t kvv64u32_get(const kvv64u32_t *v, uint64_t key, uint32_t dflt) {
  size_t L = 0, R = v->n;
  while (L < R) {
    size_t M = L + ((R - L) >> 1);
    uint64_t k = v->a[M].key;
    if (key == k)
      return v->a[M].val;
    else if (key < k)
      R = M;
    else
      L = M + 1;
  }
  return dflt;
}

/* ============================ Homopolymer Compression =======================
 */

static char *hpc_compress(const char *s) {
  size_t n = strlen(s);
  if (n == 0)
    return xstrdup("");
  char *out = (char *)xmalloc(n + 1);
  size_t w = 0;
  char prev = 0;
  for (size_t i = 0; i < n; i++) {
    char c = s[i];
    if (c >= 'a' && c <= 'z')
      c -= 32;
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      prev = 0;
      continue;
    } /* drop non-ACGT in HPC string */
    if (c != prev) {
      out[w++] = c;
      prev = c;
    }
  }
  out[w] = 0;
  return out;
}

/* ============================ Per-thread writers ============================
 */

typedef struct {
  FILE *tsv; /* per-thread TSV (no header) */
  /* optional per-class ID lists (if classify and --ids-out) */
  FILE *id_nuc, *id_org, *id_mito, *id_pltd;
  /* optional fastq writers (if classify and --emit-fastq) */
  gzFile fq_nuc, fq_org, fq_mito, fq_pltd;
} th_out_t;

static void th_open_outputs(th_out_t *O, const char *prefix, int tid,
                            int do_classify, int three_bins,
                            const char *ids_prefix, int emit_fastq) {
  memset(O, 0, sizeof(*O));
  char fn[4096];

  /* TSV part (no header), later merged */
  snprintf(fn, sizeof(fn), "%s.syncfilter.tsv.%03d", prefix, tid);
  O->tsv = fopen(fn, "w");
  if (!O->tsv) {
    perror(fn);
    die("open tsv part");
  }

  if (do_classify && ids_prefix) {
    snprintf(fn, sizeof(fn), "%s.nuclear.ids.%03d", ids_prefix, tid);
    O->id_nuc = fopen(fn, "w");
    if (three_bins) {
      snprintf(fn, sizeof(fn), "%s.mito.ids.%03d", ids_prefix, tid);
      O->id_mito = fopen(fn, "w");
      snprintf(fn, sizeof(fn), "%s.plastid.ids.%03d", ids_prefix, tid);
      O->id_pltd = fopen(fn, "w");
    } else {
      snprintf(fn, sizeof(fn), "%s.organelle.ids.%03d", ids_prefix, tid);
      O->id_org = fopen(fn, "w");
    }
  }
  if (do_classify && emit_fastq) {
    snprintf(fn, sizeof(fn), "%s.nuclear.%03d.fastq.gz", prefix, tid);
    O->fq_nuc = gzopen(fn, "wb");
    if (three_bins) {
      snprintf(fn, sizeof(fn), "%s.mito.%03d.fastq.gz", prefix, tid);
      O->fq_mito = gzopen(fn, "wb");
      snprintf(fn, sizeof(fn), "%s.plastid.%03d.fastq.gz", prefix, tid);
      O->fq_pltd = gzopen(fn, "wb");
    } else {
      snprintf(fn, sizeof(fn), "%s.organelle.%03d.fastq.gz", prefix, tid);
      O->fq_org = gzopen(fn, "wb");
    }
  }
}

static void th_close_outputs(th_out_t *O) {
  if (O->tsv)
    fclose(O->tsv);
  if (O->id_nuc)
    fclose(O->id_nuc);
  if (O->id_org)
    fclose(O->id_org);
  if (O->id_mito)
    fclose(O->id_mito);
  if (O->id_pltd)
    fclose(O->id_pltd);
  if (O->fq_nuc)
    gzclose(O->fq_nuc);
  if (O->fq_org)
    gzclose(O->fq_org);
  if (O->fq_mito)
    gzclose(O->fq_mito);
  if (O->fq_pltd)
    gzclose(O->fq_pltd);
}

static void fq_write(gzFile f, const char *name, const char *seq,
                     const char *qual) {
  gzprintf(f, "@%s\n%s\n+\n", name, seq);
  if (qual)
    gzprintf(f, "%s\n", qual);
  else {
    size_t L = strlen(seq);
    for (size_t i = 0; i < L; i++)
      gzputc(f, 'I');
    gzputc(f, '\n');
  }
}

/* =============================== Demux phase ================================
 */

static void demux_reads_round_robin(char **files, int n_files, int T,
                                    const char *prefix, char **shard_paths) {
  /* open shard gz writers */
  gzFile *out = (gzFile *)xmalloc(sizeof(gzFile) * T);
  char fn[4096];
  for (int t = 0; t < T; t++) {
    snprintf(fn, sizeof(fn), "%s.demux.%03d.fq.gz", prefix, t);
    shard_paths[t] = xstrdup(fn);
    out[t] = gzopen(fn, "wb9");
    if (!out[t]) {
      perror(fn);
      die("open demux shard");
    }
  }

  int rr = 0;
  for (int fi = 0; fi < n_files; fi++) {
    gzFile fp = gzopen(files[fi], "rb");
    if (!fp) {
      perror(files[fi]);
      die("open input reads");
    }
    kseq_t *ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
      gzprintf(out[rr], "@%s\n%s\n+\n", ks->name.s, ks->seq.s);
      if (ks->qual.l)
        gzprintf(out[rr], "%s\n", ks->qual.s);
      else {
        for (size_t i = 0; i < ks->seq.l; i++)
          gzputc(out[rr], 'I');
        gzputc(out[rr], '\n');
      }
      rr++;
      if (rr == T)
        rr = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  for (int t = 0; t < T; t++)
    gzclose(out[t]);
  free(out);
}

/* ================================ Workers ===================================
 */

typedef struct {
  int tid;
  const sm_param_t *P;
  const folds_t *F;
  int do_classify;
  int three_bins;
  int use_hpc; /* HPC for syncmer selection */
  const char *demux_path;
  th_out_t *O;
  /* quickview */
  const kvv64u32_t *mult;
  /* graph */
  const kvv64d_t *covmap;
} worker_arg_t;

/* one worker on its demux file */
static void *worker_run(void *pp) {
  worker_arg_t *A = (worker_arg_t *)pp;
  gzFile fp = gzopen(A->demux_path, "rb");
  if (!fp) {
    perror(A->demux_path);
    die("open demux shard for read");
  }
  kseq_t *ks = kseq_init(fp);

  while (kseq_read(ks) >= 0) {
    /* compute per-read proxy vector using (optional) HPC sequence */
    u64v_t sm;
    u64v_init(&sm);
    const char *seq_sync = ks->seq.s;
    char *hpc = NULL;
    int len_sync = (int)ks->seq.l;
    if (A->use_hpc) {
      hpc = hpc_compress(ks->seq.s);
      seq_sync = hpc;
      len_sync = (int)strlen(hpc);
    }
    if (len_sync >= A->P->k)
      collect_closed_syncmer_s(seq_sync, len_sync, A->P, &sm);

    double mean = 0, med = 0, q1 = 0, q3 = 0;
    int cls = 0;
    if (sm.n) {
      if (A->mult) { /* quickview: multiplicity median */
        uint32_t *cnt = (uint32_t *)xmalloc(sm.n * sizeof(uint32_t));
        for (size_t j = 0; j < sm.n; j++)
          cnt[j] = kvv64u32_get(A->mult, sm.a[j], 1u);
        qsort(cnt, sm.n, sizeof(uint32_t), cmp_u32);
        q1 = cnt[sm.n / 4];
        q3 = cnt[(3 * sm.n) / 4];
        med = (sm.n & 1) ? cnt[sm.n / 2]
                         : 0.5 * (cnt[sm.n / 2 - 1] + cnt[sm.n / 2]);
        mean = dmean_u32(cnt, sm.n);
        free(cnt);
      } else { /* graph modes: coverage median */
        double *cov = (double *)xmalloc(sm.n * sizeof(double));
        size_t m = 0;
        for (size_t j = 0; j < sm.n; j++) {
          double v;
          if (kvv64d_get(A->covmap, sm.a[j], &v))
            cov[m++] = v;
        }
        if (m) {
          uint32_t *tmp = (uint32_t *)xmalloc(m * sizeof(uint32_t));
          for (size_t j = 0; j < m; j++) {
            double v = cov[j];
            if (v < 0)
              v = 0;
            if (v > 1e9)
              v = 1e9;
            tmp[j] = (uint32_t)llround(v * 100.0);
          }
          qsort(tmp, m, sizeof(uint32_t), cmp_u32);
          q1 = tmp[m / 4] / 100.0;
          q3 = tmp[(3 * m) / 4] / 100.0;
          med = (m & 1) ? tmp[m / 2] / 100.0
                        : 0.5 * (tmp[m / 2 - 1] + tmp[m / 2]) / 100.0;
          double s = 0;
          for (size_t j = 0; j < m; j++)
            s += (double)tmp[j];
          mean = (s / (double)m) / 100.0;
          free(tmp);
        }
        free(cov);
      }
    }
    const char *cstr = "NA";
    if (A->do_classify) {
      cls = classify(med, A->F);
      cstr = class_str(cls, A->three_bins);
    }

    /* TSV part (no header) */
    fprintf(A->O->tsv, "%s\t%zu\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n", ks->name.s,
            sm.n, mean, med, q1, q3, cstr);

    /* IDs-out (optional) */
    if (A->do_classify) {
      if (A->three_bins) {
        if (cls == 0 && A->O->id_nuc)
          fprintf(A->O->id_nuc, "%s\n", ks->name.s);
        else if (cls == 1 && A->O->id_mito)
          fprintf(A->O->id_mito, "%s\n", ks->name.s);
        else if (cls == 2 && A->O->id_pltd)
          fprintf(A->O->id_pltd, "%s\n", ks->name.s);
      } else {
        if (cls == 0 && A->O->id_nuc)
          fprintf(A->O->id_nuc, "%s\n", ks->name.s);
        else if (cls != 0 && A->O->id_org)
          fprintf(A->O->id_org, "%s\n", ks->name.s);
      }
    }

    /* FASTQ bins (optional, sharded per thread) */
    if (A->do_classify) {
      if (A->three_bins) {
        if (cls == 0 && A->O->fq_nuc)
          fq_write(A->O->fq_nuc, ks->name.s, ks->seq.s,
                   (ks->qual.l ? ks->qual.s : NULL));
        else if (cls == 1 && A->O->fq_mito)
          fq_write(A->O->fq_mito, ks->name.s, ks->seq.s,
                   (ks->qual.l ? ks->qual.s : NULL));
        else if (cls == 2 && A->O->fq_pltd)
          fq_write(A->O->fq_pltd, ks->name.s, ks->seq.s,
                   (ks->qual.l ? ks->qual.s : NULL));
      } else {
        if (cls == 0 && A->O->fq_nuc)
          fq_write(A->O->fq_nuc, ks->name.s, ks->seq.s,
                   (ks->qual.l ? ks->qual.s : NULL));
        else if (cls != 0 && A->O->fq_org)
          fq_write(A->O->fq_org, ks->name.s, ks->seq.s,
                   (ks->qual.l ? ks->qual.s : NULL));
      }
    }

    if (hpc)
      free(hpc);
    u64v_free(&sm);
  }
  kseq_destroy(ks);
  gzclose(fp);
  return NULL;
}

/* ================================ Merging ===================================
 */

static void cat_file_to(FILE *dst, const char *src_path) {
  FILE *src = fopen(src_path, "r");
  if (!src) { /* silent skip for missing */
    return;
  }
  char buf[1 << 15];
  size_t n;
  while ((n = fread(buf, 1, sizeof(buf), src)) > 0) {
    if (fwrite(buf, 1, n, dst) != n) {
      perror("merge write");
      exit(2);
    }
  }
  fclose(src);
}

static void merge_tsv_parts(const char *prefix, int T) {
  char out_fn[4096];
  snprintf(out_fn, sizeof(out_fn), "%s.syncfilter.tsv", prefix);
  FILE *out = fopen(out_fn, "w");
  if (!out) {
    perror(out_fn);
    die("merge tsv");
  }
  /* header */
  fprintf(out, "read_id\tn_syncmer\tmean\tmedian\tq25\tq75\tclass\n");
  /* parts */
  char part[4096];
  for (int t = 0; t < T; t++) {
    snprintf(part, sizeof(part), "%s.syncfilter.tsv.%03d", prefix, t);
    cat_file_to(out, part);
  }
  fclose(out);
}

static void merge_ids_parts(const char *ids_prefix, int three_bins, int T) {
  const char *names2[] = {"nuclear", "organelle"};
  const char *names3[] = {"nuclear", "mito", "plastid"};
  char out_fn[4096], part[4096];
  if (three_bins) {
    for (int i = 0; i < 3; i++) {
      snprintf(out_fn, sizeof(out_fn), "%s.%s.ids", ids_prefix, names3[i]);
      FILE *out = fopen(out_fn, "w");
      if (!out) {
        perror(out_fn);
        continue;
      }
      for (int t = 0; t < T; t++) {
        snprintf(part, sizeof(part), "%s.%s.ids.%03d", ids_prefix, names3[i],
                 t);
        cat_file_to(out, part);
      }
      fclose(out);
    }
  } else {
    for (int i = 0; i < 2; i++) {
      snprintf(out_fn, sizeof(out_fn), "%s.%s.ids", ids_prefix, names2[i]);
      FILE *out = fopen(out_fn, "w");
      if (!out) {
        perror(out_fn);
        continue;
      }
      for (int t = 0; t < T; t++) {
        snprintf(part, sizeof(part), "%s.%s.ids.%03d", ids_prefix, names2[i],
                 t);
        cat_file_to(out, part);
      }
      fclose(out);
    }
  }
}

/* ================================ CLI / Main ================================
 */

static void usage(FILE *fp) {
  fprintf(
      fp,
      "Usage: syncfilter --mode {quickview|graph-early|graph-final} [options] "
      "reads.fq[.gz] ...\n\n"
      "Core:\n"
      "  --mode STR           quickview | graph-early | graph-final\n"
      "  --gfa FILE           unitig GFA (required for graph-early/final)\n"
      "  -k INT               k-mer length for closed syncmers [121]\n"
      "  -s INT               s-mer length (<=31) [27]\n"
      "  -o PREFIX            output prefix [syncfilter]\n"
      "  -t INT               worker threads [1]\n"
      "  --hpc                compute syncmers on homopolymer-compressed copy "
      "of each read\n\n"
      "Classification (optional):\n"
      "  --bin-fold A[,B]     two-bin: A splits nuc/org; three-bin: A,B split "
      "nuc/mito/plastid\n"
      "  --three-bins         enable 3-bin mode (requires two cuts for "
      "mito/plastid)\n"
      "  --ids-out PREFIX     write per-class read-id lists for post hoc "
      "filtering\n"
      "  --emit-fastq         write gz FASTQ bins (kept sharded per thread)\n\n"
      "Other:\n"
      "  --hist               write PREFIX.hist.tsv (placeholder)\n"
      "  --version, -V        print version and exit\n"
      "  --help, -h           show this help and exit\n");
}

int main(int argc, char **argv) {
  ketopt_t o = KETOPT_INIT;
  int c;
  int k = 121, s = 27, threads = 1, use_hpc = 0;
  char *mode = NULL, *gfa = NULL, *prefix = NULL, *folds_str = NULL,
       *ids_out = NULL;
  int three_bins = 0, emit_fastq = 0, write_hist = 0;

  static ko_longopt_t L[] = {{"mode", ko_required_argument, 1},
                             {"gfa", ko_required_argument, 2},
                             {"bin-fold", ko_required_argument, 3},
                             {"emit-fastq", ko_no_argument, 4},
                             {"three-bins", ko_no_argument, 5},
                             {"hist", ko_no_argument, 6},
                             {"threads", ko_required_argument, 9},
                             {"out", ko_required_argument, 10},
                             {"version", ko_no_argument, 11},
                             {"ids-out", ko_required_argument, 12},
                             {"hpc", ko_no_argument, 13},
                             {"help", ko_no_argument, 'h'},
                             {0, 0, 0}};

  while ((c = ketopt(&o, argc, argv, 1, "hk:s:o:t:V", L)) >= 0) {
    if (c == 1)
      mode = xstrdup(o.arg);
    else if (c == 2)
      gfa = xstrdup(o.arg);
    else if (c == 3)
      folds_str = xstrdup(o.arg);
    else if (c == 4)
      emit_fastq = 1;
    else if (c == 5)
      three_bins = 1;
    else if (c == 6)
      write_hist = 1;
    else if (c == 9)
      threads = atoi(o.arg);
    else if (c == 10)
      prefix = xstrdup(o.arg);
    else if (c == 11) {
      puts(SYNCFILTER_VERSION);
      return 0;
    } else if (c == 12)
      ids_out = xstrdup(o.arg);
    else if (c == 13)
      use_hpc = 1;
    else if (c == 'k')
      k = atoi(o.arg);
    else if (c == 's')
      s = atoi(o.arg);
    else if (c == 'o')
      prefix = xstrdup(o.arg);
    else if (c == 't')
      threads = atoi(o.arg);
    else if (c == 'V') {
      puts(SYNCFILTER_VERSION);
      return 0;
    } else if (c == 'h') {
      usage(stdout);
      return 0;
    } else if (c == '?') {
      fprintf(stderr, "Unknown option: %s\n", argv[o.i - 1]);
      return 1;
    } else if (c == ':') {
      fprintf(stderr, "Missing arg for: %s\n", argv[o.i - 1]);
      return 1;
    }
  }

  if (!prefix)
    prefix = xstrdup("syncfilter");
  if (!mode) {
    usage(stderr);
    return 1;
  }
  if (k <= s) {
    fprintf(stderr, "[E] require k > s\n");
    return 1;
  }
  if (o.ind >= argc) {
    fprintf(stderr, "[E] no input reads\n");
    return 1;
  }

  /* classification is optional; only if --bin-fold provided */
  folds_t F = {.f1 = 5.0, .f2 = 100.0, .three = 0};
  int do_classify = 0;
  if (folds_str) {
    char *p = folds_str;
    F.f1 = strtod(p, &p);
    if (*p == ',') {
      F.f2 = strtod(p + 1, &p);
      F.three = 1;
    }
    do_classify = 1;
  }
  if (three_bins && !F.three) {
    fprintf(stderr, "[W] --three-bins ignored (only one cut provided)\n");
    three_bins = 0;
  }
  if (!three_bins && F.three)
    three_bins = 1;

  /* ensure output directory exists (if prefix has path) */
  {
    char *dup = xstrdup(prefix);
    char *slash = strrchr(dup, '/');
    if (slash) {
      *slash = 0;
      makedir_p(dup);
    }
    free(dup);
  }
  if (ids_out) {
    char *dup = xstrdup(ids_out);
    char *slash = strrchr(dup, '/');
    if (slash) {
      *slash = 0;
      makedir_p(dup);
    }
    free(dup);
  }

  if (write_hist) {
    char hf[4096];
    snprintf(hf, sizeof(hf), "%s.hist.tsv", prefix);
    FILE *hist = fopen(hf, "w");
    if (hist) {
      fprintf(hist, "bin_center\tcount\n");
      fclose(hist);
    }
  }

  sm_param_t P = {.k = k, .s = s};

  /* Build per-mode aux data & demux first pass */
  kvv64u32_t mult;
  kvv64u32_t *Qmult = NULL;
  kvv64d_t covmap;
  kvv64d_t *Hcov = NULL;

  int T = threads > 0 ? threads : 1;
  char **shard_paths = (char **)xmalloc(sizeof(char *) * T);

  if (strcmp(mode, "quickview") == 0) {
    kvv64u32_init(&mult);
    Qmult = &mult;

    /* One pass to build multiplicity and demux simultaneously */
    /* We need multiplicity on the same selection space (HPC or not)? */
    /* For simplicity and speed, build multiplicity on original sequence (robust
       enough). If you prefer HPC-based multiplicity, demux first then a
       pre-pass over shards. */
    /* Here: build multiplicity during demux is skipped to keep demux fast; do
     * separate pass. */

    /* Separate pass to build multiplicity on original reads */
    for (int fi = o.ind; fi < argc; fi++) {
      gzFile fp = gzopen(argv[fi], "rb");
      if (!fp) {
        perror(argv[fi]);
        die("open reads");
      }
      kseq_t *ks = kseq_init(fp);
      while (kseq_read(ks) >= 0) {
        u64v_t v;
        u64v_init(&v);
        /* multiplicity is typically OK without HPC; switch below if you want
         * HPC here as well */
        collect_closed_syncmer_s(ks->seq.s, (int)ks->seq.l, &P, &v);
        if (v.n) {
          qsort(v.a, v.n, sizeof(uint64_t), cmp_u64);
          size_t m = 0;
          for (size_t i = 0; i < v.n; i++)
            if (i == 0 || v.a[i] != v.a[i - 1])
              v.a[m++] = v.a[i];
          for (size_t i = 0; i < m; i++)
            kvv64u32_push(Qmult, v.a[i], 1u);
        }
        u64v_free(&v);
      }
      kseq_destroy(ks);
      gzclose(fp);
    }
    kvv64u32_sort_merge_counts(Qmult);

    /* Now demux for threaded second pass (classification or metrics) */
    demux_reads_round_robin(argv + o.ind, argc - o.ind, T, prefix, shard_paths);

  } else if (strcmp(mode, "graph-early") == 0 ||
             strcmp(mode, "graph-final") == 0) {
    if (!gfa) {
      fprintf(stderr, "[E] --gfa required for %s\n", mode);
      return 1;
    }
    covmap = build_sm_cov_from_gfa(gfa, &P);
    Hcov = &covmap;
    demux_reads_round_robin(argv + o.ind, argc - o.ind, T, prefix, shard_paths);
  } else {
    fprintf(stderr, "[E] unknown --mode %s\n", mode);
    return 1;
  }

  /* Spawn T workers, each reads its shard and writes its own outputs */
  pthread_t *tid = (pthread_t *)xmalloc(sizeof(pthread_t) * T);
  worker_arg_t *args = (worker_arg_t *)xmalloc(sizeof(worker_arg_t) * T);
  th_out_t *outs = (th_out_t *)xmalloc(sizeof(th_out_t) * T);

  for (int t = 0; t < T; t++) {
    th_open_outputs(&outs[t], prefix, t, do_classify, three_bins, ids_out,
                    (emit_fastq && do_classify));
    args[t].tid = t;
    args[t].P = &P;
    args[t].F = &F;
    args[t].do_classify = do_classify;
    args[t].three_bins = three_bins;
    args[t].use_hpc = use_hpc;
    args[t].demux_path = shard_paths[t];
    args[t].O = &outs[t];
    args[t].mult = Qmult;
    args[t].covmap = Hcov;
    pthread_create(&tid[t], NULL, worker_run, &args[t]);
  }

  for (int t = 0; t < T; t++)
    pthread_join(tid[t], NULL);
  for (int t = 0; t < T; t++)
    th_close_outputs(&outs[t]);

  /* merge TSV parts into single headered TSV */
  merge_tsv_parts(prefix, T);

  /* merge IDs (if any) into single lists */
  if (ids_out && do_classify)
    merge_ids_parts(ids_out, three_bins, T);

  /* cleanup temp files */
  for (int t = 0; t < T; t++) {
    char part[4096];
    snprintf(part, sizeof(part), "%s.syncfilter.tsv.%03d", prefix, t);
    remove(part);
    snprintf(part, sizeof(part), "%s.demux.%03d.fq.gz", prefix, t);
    remove(part);
    if (ids_out && do_classify) {
      snprintf(part, sizeof(part), "%s.nuclear.ids.%03d", ids_out, t);
      remove(part);
      if (three_bins) {
        snprintf(part, sizeof(part), "%s.mito.ids.%03d", ids_out, t);
        remove(part);
        snprintf(part, sizeof(part), "%s.plastid.ids.%03d", ids_out, t);
        remove(part);
      } else {
        snprintf(part, sizeof(part), "%s.organelle.ids.%03d", ids_out, t);
        remove(part);
      }
    }
    free(shard_paths[t]);
  }
  free(shard_paths);
  free(outs);
  free(args);
  free(tid);

  if (Qmult) {
    free(Qmult->a);
  }
  if (Hcov) {
    free(Hcov->a);
  }

  free(mode);
  free(gfa);
  free(prefix);
  free(folds_str);
  free(ids_out);
  return 0;
}

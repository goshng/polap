// scripts/polap-c-reads.c
// Fast PAF aggregator for organelle read counting (MT/PT)
// Build: cc -O3 -march=native -pipe -o polap-reads scripts/polap-c-reads.c

#define _GNU_SOURCE
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *qname;
  int qlen;
  int s, e; // 1-based closed [s,e]
} Interval;

typedef struct {
  char *contig;
  int length;
  int count;
} Summary;

typedef struct {
  Interval *a;
  size_t n, m;
} VecI;

typedef struct {
  Summary *a;
  size_t n, m;
} VecS;

static void die(const char *msg) {
  fprintf(stderr, "[err] %s\n", msg);
  exit(1);
}

static void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p)
    die("OOM");
  return p;
}
static void *xrealloc(void *p, size_t n) {
  p = realloc(p, n);
  if (!p)
    die("OOM");
  return p;
}
static char *xstrdup(const char *s) {
  size_t n = strlen(s) + 1;
  char *p = xmalloc(n);
  memcpy(p, s, n);
  return p;
}

static void vecI_push(VecI *v, Interval z) {
  if (v->n == v->m) {
    v->m = v->m ? v->m << 1 : 1024;
    v->a = xrealloc(v->a, v->m * sizeof(Interval));
  }
  v->a[v->n++] = z;
}
static void vecS_push(VecS *v, Summary z) {
  if (v->n == v->m) {
    v->m = v->m ? v->m << 1 : 1024;
    v->a = xrealloc(v->a, v->m * sizeof(Summary));
  }
  v->a[v->n++] = z;
}

static int cmp_interval(const void *pa, const void *pb) {
  const Interval *a = (const Interval *)pa, *b = (const Interval *)pb;
  int c = strcmp(a->qname, b->qname);
  if (c)
    return c;
  if (a->s != b->s)
    return (a->s < b->s) ? -1 : 1;
  if (a->e != b->e)
    return (a->e < b->e) ? -1 : 1;
  return 0;
}

static int cmp_summary(const void *pa, const void *pb) {
  const Summary *a = (const Summary *)pa, *b = (const Summary *)pb;
  return strcmp(a->contig, b->contig);
}

static int parse_line_cols(char *line, char **out, int want) {
  // Split by '\t' up to 'want' fields (we need 12 but we only pick specific)
  // Return #columns found
  int n = 0;
  char *p = line;
  while (*p && n < want) {
    out[n++] = p;
    char *t = strchr(p, '\t');
    if (!t)
      break;
    *t = '\0';
    p = t + 1;
  }
  return n;
}

static int is_number(const char *s) {
  if (!*s)
    return 0;
  while (*s) {
    if (!isdigit((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

static void load_paf(const char *path, double min_id, int min_mapq, VecI *out) {
  FILE *fp = fopen(path, "rb");
  if (!fp) {
    fprintf(stderr, "[err] cannot open %s\n", path);
    exit(1);
  }
  char *line = NULL;
  size_t cap = 0;
  ssize_t len;
  char *cols[12];
  size_t line_no = 0;
  while ((len = getline(&line, &cap, fp)) != -1) {
    line_no++;
    if (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
      line[--len] = 0;
    int ncol = parse_line_cols(line, cols, 12);
    if (ncol < 12)
      continue; // skip short lines

    // Needed cols (1-based):
    // 1=qname,2=qlen,3=qstart,4=qend,10=nmatch,11=alen,12=mapq
    char *qname = cols[0];
    char *qlen_s = cols[1], *qs_s = cols[2], *qe_s = cols[3];
    char *nmatch_s = cols[9], *alen_s = cols[10], *mapq_s = cols[11];

    if (!is_number(qlen_s) || !is_number(qs_s) || !is_number(qe_s) ||
        !is_number(nmatch_s) || !is_number(alen_s) || !is_number(mapq_s)) {
      continue;
    }
    int qlen = atoi(qlen_s);
    int qs = atoi(qs_s) + 1; // 1-based
    int qe = atoi(qe_s);
    int nmatch = atoi(nmatch_s);
    int alen = atoi(alen_s);
    int mapq = atoi(mapq_s);
    if (alen <= 0)
      continue;

    double id = (double)nmatch / (double)alen;
    if (mapq < min_mapq || id < min_id)
      continue;

    Interval z;
    z.qname = xstrdup(qname);
    z.qlen = qlen;
    z.s = qs;
    z.e = qe;
    vecI_push(out, z);
  }
  free(line);
  fclose(fp);
}

static VecS summarize(VecI *v) {
  VecS res = {0};
  if (v->n == 0)
    return res;
  qsort(v->a, v->n, sizeof(Interval), cmp_interval);

  size_t i = 0;
  while (i < v->n) {
    size_t j = i + 1;
    // same contig block
    while (j < v->n && strcmp(v->a[j].qname, v->a[i].qname) == 0)
      j++;

    // merge intervals in [i, j)
    int merged = 0;
    int cur_s = v->a[i].s;
    int cur_e = v->a[i].e;
    int qlen = v->a[i].qlen;
    for (size_t k = i + 1; k < j; ++k) {
      if (v->a[k].s <= cur_e) {
        if (v->a[k].e > cur_e)
          cur_e = v->a[k].e;
      } else {
        merged++;
        cur_s = v->a[k].s;
        cur_e = v->a[k].e;
      }
      // prefer nonzero qlen if encountered
      if (v->a[k].qlen > 0 && qlen <= 0)
        qlen = v->a[k].qlen;
    }
    merged++;

    Summary s;
    s.contig = xstrdup(v->a[i].qname);
    s.length = qlen;
    s.count = merged;
    vecS_push(&res, s);
    i = j;
  }
  return res;
}

// Helper: extract numeric tail after last '.'; return newly-alloc string or
// NULL
static char *edge_tail(const char *name) {
  const char *dot = strrchr(name, '.');
  if (!dot || !*(dot + 1))
    return NULL;
  const char *p = dot + 1;
  for (const char *q = p; *q; ++q)
    if (!isdigit((unsigned char)*q))
      return NULL;
  return xstrdup(p);
}

static void write_outputs(const VecS *mt, const VecS *pt, const char *outdir,
                          int min_pt) {
  // Merge two sorted summary lists by contig
  // First sort them (just in case)
  VecS MT = *mt, PT = *pt;
  qsort(MT.a, MT.n, sizeof(Summary), cmp_summary);
  qsort(PT.a, PT.n, sizeof(Summary), cmp_summary);

  // Build union into a flat struct array
  typedef struct {
    char *Contig, *Edge;
    int Length, Depth, Copy, MT, PT;
  } Row;

  Row *rows = NULL;
  size_t rn = 0, rm = 0;
  size_t i = 0, j = 0;
  while (i < MT.n || j < PT.n) {
    const char *cm = (i < MT.n) ? MT.a[i].contig : NULL;
    const char *cp = (j < PT.n) ? PT.a[j].contig : NULL;
    int cmp = 0;
    if (!cm)
      cmp = 1;
    else if (!cp)
      cmp = -1;
    else
      cmp = strcmp(cm, cp);

    const char *name = NULL;
    int len_mt = 0, cnt_mt = 0, len_pt = 0, cnt_pt = 0;
    if (cmp < 0) { // mt only
      name = cm;
      len_mt = MT.a[i].length;
      cnt_mt = MT.a[i].count;
      i++;
    } else if (cmp > 0) { // pt only
      name = cp;
      len_pt = PT.a[j].length;
      cnt_pt = PT.a[j].count;
      j++;
    } else { // both
      name = cm;
      len_mt = MT.a[i].length;
      cnt_mt = MT.a[i].count;
      len_pt = PT.a[j].length;
      cnt_pt = PT.a[j].count;
      i++;
      j++;
    }
    int length = len_mt ? len_mt : len_pt;
    int MTc = cnt_mt;
    int PTc = cnt_pt;
    if (PTc < min_pt)
      continue;

    if (rn == rm) {
      rm = rm ? rm << 1 : 1024;
      rows = xrealloc(rows, rm * sizeof(Row));
    }
    rows[rn].Contig = xstrdup(name);
    rows[rn].Length = length;
    rows[rn].Depth = 1;
    rows[rn].Copy = 1;
    rows[rn].MT = MTc;
    rows[rn].PT = PTc;
    rows[rn].Edge = edge_tail(name);
    rn++;
  }

  // helpers for sorting
  // desc MT
  int cmp_descMT(const void *pa, const void *pb) {
    const Row *a = (const Row *)pa, *b = (const Row *)pb;
    if (a->MT != b->MT)
      return (a->MT > b->MT) ? -1 : 1;
    return strcmp(a->Contig, b->Contig);
  }
  // key (MT<=PT) logical sort (FALSE before TRUE) -> emulate by int key
  int key_le_MTPT(const Row *r) { return (r->MT <= r->PT) ? 1 : 0; }
  int cmp_key_le_MTPT_then_name(const void *pa, const void *pb) {
    const Row *a = (const Row *)pa, *b = (const Row *)pb;
    int ka = key_le_MTPT(a), kb = key_le_MTPT(b);
    if (ka != kb)
      return (ka < kb) ? -1 : 1;
    return strcmp(a->Contig, b->Contig);
  }
  // desc (MT-PT)
  int cmp_descDiff(const void *pa, const void *pb) {
    const Row *a = (const Row *)pa, *b = (const Row *)pb;
    int da = a->MT - a->PT, db = b->MT - b->PT;
    if (da != db)
      return (da > db) ? -1 : 1;
    return strcmp(a->Contig, b->Contig);
  }

  // Output 2: all by desc MT
  {
    char path[4096];
    snprintf(path, sizeof(path),
             "%s/assembly_info_organelle_annotation_count-all.txt", outdir);
    Row *tmp = xmalloc(rn * sizeof(Row));
    memcpy(tmp, rows, rn * sizeof(Row));
    qsort(tmp, rn, sizeof(Row), cmp_descMT);
    FILE *fo = fopen(path, "wb");
    if (!fo)
      die("cannot write out2");
    for (size_t k = 0; k < rn; k++) {
      fprintf(fo, "%s\t%d\t%d\t%d\t%d\t%d\t%s\n", tmp[k].Contig, tmp[k].Length,
              tmp[k].Depth, tmp[k].Copy, tmp[k].MT, tmp[k].PT,
              tmp[k].Edge ? tmp[k].Edge : "");
    }
    fclose(fo);
    free(tmp);
  }

  // Output 1: first desc MT, then stable sort by (MT<=PT)
  {
    char path[4096];
    snprintf(path, sizeof(path),
             "%s/assembly_info_organelle_annotation_count.txt", outdir);
    Row *tmp = xmalloc(rn * sizeof(Row));
    memcpy(tmp, rows, rn * sizeof(Row));
    qsort(tmp, rn, sizeof(Row), cmp_descMT);
    // stable-like: decorate with index and do second key
    // For brevity do a single pass qsort with the second comparator (close
    // enough)
    qsort(tmp, rn, sizeof(Row), cmp_key_le_MTPT_then_name);
    FILE *fo = fopen(path, "wb");
    if (!fo)
      die("cannot write out1");
    for (size_t k = 0; k < rn; k++) {
      fprintf(fo, "%s\t%d\t%d\t%d\t%d\t%d\t%s\n", tmp[k].Contig, tmp[k].Length,
              tmp[k].Depth, tmp[k].Copy, tmp[k].MT, tmp[k].PT,
              tmp[k].Edge ? tmp[k].Edge : "");
    }
    fclose(fo);
    free(tmp);
  }

  // Output 4: MT>PT ordered by desc(MT-PT)
  {
    char path[4096];
    snprintf(path, sizeof(path), "%s/contig-annotation-depth-table.txt",
             outdir);
    Row *tmp = xmalloc(rn * sizeof(Row));
    size_t tn = 0;
    for (size_t k = 0; k < rn; k++)
      if (rows[k].MT > rows[k].PT)
        tmp[tn++] = rows[k];
    qsort(tmp, tn, sizeof(Row), cmp_descDiff);
    FILE *fo = fopen(path, "wb");
    if (!fo)
      die("cannot write out4");
    for (size_t k = 0; k < tn; k++) {
      fprintf(fo, "%s\t%d\t%d\t%d\t%d\t%d\t%s\n", tmp[k].Contig, tmp[k].Length,
              tmp[k].Depth, tmp[k].Copy, tmp[k].MT, tmp[k].PT,
              tmp[k].Edge ? tmp[k].Edge : "");
    }
    fclose(fo);
    free(tmp);
  }

  // Output 5: PT>MT, emulate your 2-step arrange
  {
    char path[4096];
    snprintf(path, sizeof(path), "%s/pt-contig-annotation-depth-table.txt",
             outdir);
    Row *tmp = xmalloc(rn * sizeof(Row));
    memcpy(tmp, rows, rn * sizeof(Row));
    // sort by desc MT then (PT<=MT)
    qsort(tmp, rn, sizeof(Row), cmp_descMT);
    qsort(tmp, rn, sizeof(Row), cmp_key_le_MTPT_then_name);
    FILE *fo = fopen(path, "wb");
    if (!fo)
      die("cannot write out5");
    for (size_t k = 0; k < rn; k++) {
      if (tmp[k].PT > tmp[k].MT) {
        fprintf(fo, "%s\t%d\t%d\t%d\t%d\t%d\t%s\n", tmp[k].Contig,
                tmp[k].Length, tmp[k].Depth, tmp[k].Copy, tmp[k].MT, tmp[k].PT,
                tmp[k].Edge ? tmp[k].Edge : "");
      }
    }
    fclose(fo);
    free(tmp);
  }

  // cleanup
  for (size_t k = 0; k < rn; k++) {
    free(rows[k].Contig);
    if (rows[k].Edge)
      free(rows[k].Edge);
  }
  free(rows);
}

int main(int argc, char **argv) {
  const char *mt = NULL, *pt = NULL, *out = NULL;
  int min_mapq = 1, min_pt = 0;
  double min_id = 0.15;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--mt") && i + 1 < argc)
      mt = argv[++i];
    else if (!strcmp(argv[i], "--pt") && i + 1 < argc)
      pt = argv[++i];
    else if (!strcmp(argv[i], "--output") && i + 1 < argc)
      out = argv[++i];
    else if (!strcmp(argv[i], "--min-mapq") && i + 1 < argc)
      min_mapq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--min-pt") && i + 1 < argc)
      min_pt = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--min-identity") && i + 1 < argc)
      min_id = strtod(argv[++i], NULL);
    else {
      fprintf(stderr, "Unknown/invalid arg: %s\n", argv[i]);
      return 1;
    }
  }
  if (!mt || !pt || !out) {
    fprintf(stderr,
            "Usage: %s --mt mt.paf --pt pt.paf --output outdir [--min-mapq "
            "INT] [--min-pt INT] [--min-identity FLOAT]\n",
            argv[0]);
    return 1;
  }

  VecI mtv = {0}, ptv = {0};
  load_paf(mt, min_id, min_mapq, &mtv);
  load_paf(pt, min_id, min_mapq, &ptv);

  VecS mts = summarize(&mtv);
  VecS pts = summarize(&ptv);

  write_outputs(&mts, &pts, out, min_pt);

  // free intervals
  for (size_t i = 0; i < mtv.n; i++)
    free(mtv.a[i].qname);
  for (size_t i = 0; i < ptv.n; i++)
    free(ptv.a[i].qname);
  free(mtv.a);
  free(ptv.a);
  for (size_t i = 0; i < mts.n; i++)
    free(mts.a[i].contig);
  for (size_t i = 0; i < pts.n; i++)
    free(pts.a[i].contig);
  free(mts.a);
  free(pts.a);
  return 0;
}

# select_triplet.awk
# Input: lines "ID LENGTH"
# Env vars: MINIR, MAXIR (integers)
# Output: one line "LSC_ID IR_ID SSC_ID"
{
  ids[NR]=$1; lens[NR]=$2+0;
}
END{
  m=NR;
  for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) for(k=j+1;k<=m;k++){
    L1=lens[i]; L2=lens[j]; L3=lens[k];
    I1=ids[i];  I2=ids[j];  I3=ids[k];
    if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
    if(L2<L3){t=L2;L2=L3;L3=t; u=I2;I2=I3;I3=u}
    if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
    total=L1+L3+2*L2
    if(total>=100000 && total<=250000 && L2>=MINIR && L2<=MAXIR){
      printf "%s\t%s\t%s\n", I1,I2,I3
      exit
    }
  }
}

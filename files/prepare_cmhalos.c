#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {

  FILE * fin, * fout;
  char buf[300];
  int i, j, k, ninput_files, cont, dummy, nhalos, halosize=13;
  float * haloblock;

  if(argc < 3) {
    fprintf(stdout, "\nParameters are missing.\n");
    fprintf(stdout, "Call with <Input_Filename, Number_of_Input_Files>\n\n");
    exit(0);
  }

  ninput_files = atoi(argv[2]);

  for (i=0; i<ninput_files-1; i++) {

    sprintf(buf, "%s.%d", argv[1], i);
    fin  = fopen(buf, "rb");
    printf("Opened: %s\n", buf); 

    sprintf(buf, "%s_reform_c.%d", argv[1], i);
    fout = fopen(buf, "w");

    fread(&dummy, sizeof(dummy), 1, fin);
    fread(&nhalos, sizeof(int), 1, fin);
    fread(&dummy, sizeof(dummy), 1, fin);
    fprintf(fout, "%13d\n", nhalos);
    haloblock = (float*)malloc(nhalos*halosize*sizeof(float));
    fread(&dummy, sizeof(dummy), 1, fin);
    fread(&(haloblock[0]), sizeof(float), nhalos*halosize, fin);
    fread(&dummy, sizeof(dummy), 1, fin);
    for (j=0; j<nhalos; j++) {
      fprintf(fout, "%13d\t", (int)haloblock[j*halosize]);
      for (k=1; k<halosize; k++) fprintf(fout, "%13.6lf\t", haloblock[j*halosize+k]);
      fprintf(fout, "\n");
    }
    fclose(fin);
    fclose(fout);
  }

  return 0;
}

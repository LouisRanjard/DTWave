/*
 =========================================================================================================
 dtwave.c
 
 2008-2016
 Louis Ranjard
 Compute dynamic time warp distance and a weighted average vector sequence.
 
 
Copyright (C) 2016 Louis Ranjard

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 
 dtwave [-c | -w NUM | -a] [DIR | FILE FILE [FILE]]
 
 -c
 use compression/expansion
 
 -w NUM
 use weight "NUM" to compute average sequence, use 0.5 for equal weights of both input sequences
 
 -a DIR
 compute the average sequence of a whole directory rather than just the pairwise distance matrix 
 
 
 Preprocessor Directives:
 
 -DBIT64
 compile on a 64bits architecture
 
 -DMACOS
 compile on macos 64bits architecture (use with -m64)
 
 -DDEBUG
 print debugging messages
 
 
 Created by Louis Ranjard <louis.ranjard@gmail.com> on 07/08/08.
 
 
 
 Version history:
 v2016-05-10
- speed up
 v2014-05-16
 - windows executable now available
 v2014-03-20
 - error messages are now more informative
 v2013-05-19
 - portable to 64bits system
 - bug fix in the average sequence computation (averseq, warpav)
 - compression/expansion are now optionals (-c)
 
 =========================================================================================================
 */

/*
compilation: 
gcc -g -DDEBUG dtwave.c -lm (also potentialy -DBIT64, -DMACOS, -m64 or -m32)
valgrind --tool=memcheck --leak-check=yes ./a.out ./test_files/sound1.parameters ./test_files/sound2.parameters
for production:
gcc -O3 -m64 -DMACOS dtwave.c -o dtwave_mac -lm
gcc -O3 -m64 -DBIT64 dtwave.c -o dtwave_unix64 -lm
gcc -O3 -m32 dtwave.c -o dtwave_unix32 -lm
Windows (first install MinGW then): gcc -O3 dtwave.c -o dtwave_win.exe

usage:
dtwave ./test_files/sound1.parameters ./test_files/sound2.parameters
dtwave -c ./test_files/sound1.parameters ./test_files/sound2.parameters
dtwave -w .6 ./test_files/sound1.parameters ./test_files/sound2.parameters ./test_files/average6.parameters
dtwave ./test_files
dtwave -a ./test_files
*/

// need to link the math library for compiling gcc -lm dtwave.c
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <sys/stat.h> 
#include <dirent.h> 
#include <string.h>

#ifdef BIT64
typedef struct {/* macosx 64bits */
	int nSamples;
	int sampPeriod;
	short sampSize;
	short parmKind;
} htk_header_t;
#else
typedef struct {/* linux 32bits */
	long nSamples;
	long sampPeriod;
	short sampSize;
	short parmKind;
} htk_header_t;
#endif

typedef struct  {
	int numfiles;
	char **filenam;
} dir_t;

typedef struct  {
	htk_header_t *header1;
	htk_header_t *header2;
	float *mat1;
	float *mat2;
	float *distance;
	float *cow;
	float *d;
	float *traceb;
	float *indelc;
	float *mataverage;
	float weight;
	int *duree;
	int numrows;
	int numcols;
	int numv;
	char *nameave;
} align_t;

/* GLOBAL VARIABLES */
htk_header_t **ADDRESSES_HTK_HEADER_T; // store all used addresses
int LASTADD_HTK_HEADER_T=0; // store how many have actually been used
float **ADDRESSES_FLOAT; // store all used addresses
int LASTADD_FLOAT=0; // store how many have actually been used
int **ADDRESSES_INT; // store all used addresses
int LASTADD_INT=0; // store how many have actually been used
char **ADDRESSES_CHAR; // store all used addresses
int LASTADD_CHAR=0; // store how many have actually been used
int SWAP; // endianality

/* PROTOTYPES */
float warpav_ce(float *traceb, float *d, int numrows, int numcols, float *indelc);
float warpav(float *traceb, float *d, int numrows, int numcols, float *indelc);
float *averseq(int *duree, float *traceb, float *mat1, float *mat2, float weight, int numrows, int numcols, int numv);
int mincost(float cost[], int a, int b);
int maxcost(float cost[], int a, int b);
float louiround(float x);
void euclidist(float *d, float *mat1, float *mat2, float *indelc, int numrows, int numcols, int numv, float *cow);
int CheckArchEndianality(void);
void SwapLong32(long *p);
void SwapInt(int *p);
void SwapShort(short *p);
void SwapFloat(float *p);
float *ReadHTKFile(char *name, int *size1, int *size2, htk_header_t *header);
void WriteHTKFile(char *name, float *data, long nSamples, long sampPeriod, short sampSize, short parmKind);
int isDir(const char* target);
void listdir(char *name, dir_t *directo);
void UTMatIndex(int N, int idx, int *n, int *m);
int VectUIMatIndex(int N, int n, int m);
int npairs(int nitem);
void Align(align_t *alignment, dir_t directo, int id1, int id2, int ce);
void free_align_t(align_t *alignment);
void null_align_t(align_t *alignment);
void print_align_t(align_t *alignment);


float warpav_ce(float *traceb, float *d, int numrows, int numcols, float *indelc)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    float cost[5] ;
    float *dsf ;
    float dist ;
	
    dsf = (float *)malloc((numrows+1) * (numcols+1) * sizeof(float));
    assert(dsf!=NULL);
	
    /* first cell (1,1) */
    *dsf = *d * 2;
    *traceb = 0;
	
    /* first column */
	for (i=1;i<=numcols;i++) {
		*(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
		*(traceb+((numrows+1)*i)) = 1;
	}
	
    /* first row */
	for (j=1;j<=numrows;j++) {
		*(dsf+j) = *(dsf+(j-1)) + *indelc;
		*(traceb+j) = 3;
	}
	
    /* cell (2,2) */
	cost[1] = *(dsf+1) + *indelc;
	cost[2] = *dsf + *d;
	cost[3] = *(dsf+numrows+1) + *indelc;
	k = mincost(cost,1,3);
	*(traceb+numrows+2) = k;
	*(dsf+numrows+2) = cost[k];
	
    /* second column */
	j=1;
	for (i=2;i<=numcols;i++) {
		cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
		cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
		cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
		cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
		k = mincost(cost,0,3);
		*(traceb+((numrows+1)*i)+1) = k;
		*(dsf+((numrows+1)*i)+1) = cost[k];
	}
	
    /* second row */
	i=1;
	for (j=2;j<=numrows;j++) {
		cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
		cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
		cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
		cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
		k = mincost(cost,1,4);
		*(traceb+((numrows+1))+j) = k;
		*(dsf+(numrows+1)+j) = cost[k];
	}
	
    /* rest of the matrix */
	for (i=2;i<=numcols;i++) {
		for (j=2;j<=numrows;j++) {
			cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
			cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
			cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
			cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
			cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
			k = mincost(cost,0,4);
			*(traceb+((numrows+1)*i)+j) = k;
			*(dsf+((numrows+1)*i)+j) = cost[k];
		}
	}
	
	/* compute distance */
	dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (float)(numrows+numcols);
	/* end free space
	 double mini=*(dsf+((numrows+1)*(numcols+1))-1), minj=*(dsf+((numrows+1)*(numcols+1))-1) ;
	 for (i=0;i<=numcols;i++) {if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {mini=*(dsf+((numrows+1)*(i))+numrows);} }
	 for (j=0;j<=numrows;j++) {if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {minj=*(dsf+((numrows+1)*(numcols))+j);} }
	 if (minj<mini) {dist = minj/(double)(numrows+numcols);} else {dist = mini/(double)(numrows+numcols);} */
	
	free(dsf);
	return dist ;
}

float warpav(float *traceb, float *d, int numrows, int numcols, float *indelc)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    float cost[4] ; /* need to keep 4 cells for compatibility with warpav_ce */
    float *dsf ;
    float dist ;
	
    dsf = (float *)malloc((numrows+1) * (numcols+1) * sizeof(float));
    assert(dsf!=NULL);
	
    /* first cell (1,1) */
    *dsf = *d * 2;
    *traceb = 0;
	
    /* first column */
	for (i=1;i<=numcols;i++) {
		*(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
		*(traceb+((numrows+1)*i)) = 1;
	}
	
    /* first row */
	for (j=1;j<=numrows;j++) {
		*(dsf+j) = *(dsf+(j-1)) + *indelc;
		*(traceb+j) = 3;
	}
	
    /* rest of the matrix */
	for (i=1;i<=numcols;i++) {
		for (j=1;j<=numrows;j++) {
			cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
			cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
			cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
			k = mincost(cost,1,3);
			*(traceb+((numrows+1)*i)+j) = k;
			*(dsf+((numrows+1)*i)+j) = cost[k];
		}
	}
	
	/* compute distance */
	dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (float)(numrows+numcols);
	/* end free space
	 double mini=*(dsf+((numrows+1)*(numcols+1))-1), minj=*(dsf+((numrows+1)*(numcols+1))-1) ;
	 for (i=0;i<=numcols;i++) {if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {mini=*(dsf+((numrows+1)*(i))+numrows);} }
	 for (j=0;j<=numrows;j++) {if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {minj=*(dsf+((numrows+1)*(numcols))+j);} }
	 if (minj<mini) {dist = minj/(double)(numrows+numcols);} else {dist = mini/(double)(numrows+numcols);} */
	
	free(dsf);
	return dist ;
}

/* compute the average sequence given the alignment */
float *averseq(int *duree, float *traceb, float *mat1, float *mat2, float weight, int numrows, int numcols, int numv)
{
    int i,j,n=0;
    int l=0;
    int maxt=0,maxtmp=0; /* length of the matrix */
    int ta,tb;
    float *mat4,*t,*mata;
	
    mat4 = (float *)malloc((numcols+numrows)*numv*sizeof(float)) ;
    assert(mat4!=NULL);
    t = (float *)malloc((numcols+numrows)*sizeof(float)) ;
    assert(t!=NULL);
	
    i = numcols; //printf(" i=%i\n",i );
    j = numrows; //printf(" j=%i\n",j );
	
    /* trace back the alignment and implement the average vector sequence */
    while ( i>0 || j>0 ) {
		  *(t+l) = weight*j + (1-weight)*i; //printf(" %i , %i -> %f(%f) w=%f\n",i,j,*(t+l),louiround(*(t+l)),weight );
		  maxtmp = (int)(louiround(*(t+l)));
		  if ( maxtmp>maxt ) { maxt=maxtmp; }
		  switch ( (int) *(traceb+((numrows+1)*(i))+(j)) ) {
			  case 0:
				  for (n=0;n<numv;n++) {
					  *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-2)))) ; }
				  i-=2;
				  j--;
				  break;
			  case 1:
				  for (n=0;n<numv;n++) {
					  *(mat4+n+(numv*l)) = (*(mat2+n+(numv*(i-1)))); }
				  i--;
				  break;
        case 2:
           for (n=0;n<numv;n++) { 
                *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
           i--;
           j--;
           break;
			  case 3:
				  for (n=0;n<numv;n++) {
					  *(mat4+n+(numv*l)) = (*(mat1+n+(numv*(j-1)))); }
				  j--;
				  break;
			  case 4:
				  for (n=0;n<numv;n++) {
					  *(mat4+n+(numv*l)) = 0.5*weight*(*(mat1+n+(numv*(j-1)))) + 0.5*weight*(*(mat1+n+(numv*(j-2)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
				  i--;
				  j-=2;
				  break;
			  default:
				  printf("error: traceback no match\n");
          break;
		  }
		  l++;
    }
	
    mata = (float *)malloc((maxt)*numv*sizeof(float)) ;
    assert(mata!=NULL);
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = mata;

    /* compute average with correct time direction and number of time points*/
    for ( ta=1 ; ta<=maxt ; ta++ ) {
		for ( tb=l-1 ; tb>=0 ; tb-- ) {
			if ( *(t+tb)==(float)ta ) {
				for ( n=0 ; n<numv ; n++ ) {
					*(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
				break;
			} else if ( *(t+tb)>(float)ta ) {
				if ( tb==l-1 ) {
					for ( n=0 ; n<numv ; n++ ) {
						*(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
				} else {
					for ( n=0 ; n<numv ; n++ ) {
						*(mata+n+(numv*(ta-1))) = 0.5*(*(mat4+n+(numv*(tb+1)))) + 0.5*(*(mat4+n+(numv*(tb)))) ; }
				}
				break;
			} else if ( tb==0 ) {
				for ( n=0 ; n<numv ; n++ ) { 
					*(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
			}
		}
    }
	
    *duree = maxt ;
    free(mat4);
    free(t);
    return mata ;
}

/* find the minimum value in an array between indexes a and b */
int mincost(float cost[], int a, int b)
{
    /*int j=a;
	 while (a<b+1) {
	 if (cost[a]<cost[j]) {j=a;}
	 a++;
	 }
	 return j;*/
    /* choose one randomly if equals */
    int j=a;
    while (a<b+1) {
        if (cost[a]<cost[j]) {j=a;}
        else if (cost[a]==cost[j]) {
		if (time(NULL)%2==0) {j=a;} }
        a++;
    }
    return j;
}

/* find the minimum value in an array between indexes a and b */
int maxcost(float cost[], int a, int b)
{   /* choose one randomly if equals */
    int j=a;
    while (a<b+1) {
        if (cost[a]>cost[j]) {j=a;}
        else if (cost[a]==cost[j]) {
		if (time(NULL)%2==0) {j=a;} }
        a++;
    }
    return j;
}

/* no round() function in math.h so here it is, SOMETIMES THIS FUNCTION IS BUILT IN */
float louiround(float x)
{
    double intpart;
    if(modf(x,&intpart)>=0.5){
		return ceil(x);}
    else{
		return floor(x);}
}

/* compute Euclidean distance between vectors of 2 matrices 
 the insertion deletion cost is set as (half) the average distance between vectors */
void euclidist(float *d, float *mat1, float *mat2, float *indelc, int numrows, int numcols, int numv, float *cow)
{
    int i,j,n ;
	
    for (i=0;i<numcols;i++) {
		for (j=0;j<numrows;j++) {
			*(d+((numrows)*i)+j)=0;
			for ( n=0 ; n<numv ; n++ ) {
				*(d+((numrows)*i)+j) = *(d+((numrows)*i)+j) + *(cow+n) * ((*(mat2+n+(numv*(i))))-(*(mat1+n+(numv*(j)))))* ((*(mat2+n+(numv*(i))))-(*(mat1+n+(numv*(j))))) ; }
			*(d+((numrows)*i)+j) = sqrt( *(d+((numrows)*i)+j) ) ;
			*indelc += *(d+(numrows*i)+j) ;
		}
    }
    *indelc = *indelc/(numrows*numcols) ;
    return ;
}

int CheckArchEndianality()
{
    int Endian = 0x00000001; // assuming target architecture is 32-bit
    // as Endian = 0x00000001 so MSB (Most Significant Byte) = 0x00 and LSB (Least Significant Byte) = 0x01    
	// casting down to a single byte value LSB discarding higher bytes	
    return (*(char *) &Endian == 0x01) ? 1 : 0; // LITTLEENDIAN:1, BIGENDIAN:0
}

/* swap byte order of long data value *p */
void SwapLong32(long *p)
{
	char temp,*q;
	q = (char*) p;
	temp = *q; *q = *(q+3); *(q+3) = temp;
	temp = *(q+1); *(q+1) = *(q+2); *(q+2) = temp;
}

/* swap byte order of int data value *p */
void SwapInt(int *p)
{
	char temp,*q;
	q = (char*) p;
	temp = *q; *q = *(q+3); *(q+3) = temp;
	temp = *(q+1); *(q+1) = *(q+2); *(q+2) = temp;
}

/* swap byte order of short data value *p */
void SwapShort(short *p)
{
	char temp,*q;
	q = (char*) p;
	temp = *q; *q = *(q+1); *(q+1) = temp;
}

/* swap byte order of float data value *p */
void SwapFloat(float *p)
{
	char temp,*q;
	q = (char*) p;
	temp = *q; *q = *(q+3); *(q+3) = temp;
	temp = *(q+1); *(q+1) = *(q+2); *(q+2) = temp;
}

float *ReadHTKFile(char *name, int *size1, int *size2, htk_header_t *header)
{
	FILE *file;
	float *buffer;
	float *mat;
	unsigned long fileLen;
	int i,j;
	
	//Open file
	file = fopen(name, "rb");
	if (!file) {
		fprintf(stderr, "Unable to open file %s\n", name);
		return 0; }
	else { fprintf(stderr, "read htk file %s\n", name); }
	
	//Get file length
	fseek(file, 0, SEEK_END);
	fileLen = (ftell(file)-sizeof(htk_header_t))/4; //printf("file length=%i\n", fileLen);
	fseek(file, 0, SEEK_SET);
	
	//read HTK header
	j = fread(header, sizeof(htk_header_t), 1, file);
  #ifdef BIT64
	if (SWAP) {SwapInt(&(*header).nSamples);} //printf("ns=%i\n", (*header).nSamples);
	if (SWAP) {SwapInt(&(*header).sampPeriod);} //printf("sp=%i\n", (*header).sampPeriod);
  #else
	if (SWAP) {SwapLong32(&(*header).nSamples);} //printf("ns=%i\n", (*header).nSamples);
	if (SWAP) {SwapLong32(&(*header).sampPeriod);} //printf("sp=%i\n", (*header).sampPeriod);
	#endif
	if (SWAP) {SwapShort(&(*header).sampSize);} //printf("ss=%hi\n", (*header).sampSize);
	if (SWAP) {SwapShort(&(*header).parmKind);} //printf("pk=%hi\n", (*header).parmKind);
	*size1 = (int)(*header).nSamples;
	*size2 = (int)((*header).sampSize/4);
	
	//Allocate memory
	buffer=(float *)malloc(fileLen*sizeof(float));
        assert(buffer!=NULL);
	mat = (float *)malloc((*size1)*(*size2)*sizeof(float));
        assert(mat!=NULL);
	*(ADDRESSES_FLOAT + LASTADD_FLOAT++) = mat;

	if (!buffer) {
		fprintf(stderr, "Memory error!");
		fclose(file);
		return 0;}
	
	//Read file contents into buffer
	j = fread(buffer, sizeof(float), fileLen, file);
	fclose(file);
	
	//swap bytes if required (endianality?)
	for(i = 0; i < fileLen; ++i) {
		if (SWAP) { SwapFloat(&((float *)buffer)[i]); }
		*(mat+i) = *&((float *)buffer)[i];
		//printf("%f ", ((float *)buffer)[i]);
	}
	
	free(buffer);
	return mat;
}

void WriteHTKFile(char *name, float *data, long nSamples, long sampPeriod, short sampSize, short parmKind)
{
	FILE *file;
	htk_header_t header;
	unsigned long fileLen, totalparam;
	int i;
	float *buffer;
	
	//Open file
	file = fopen(name, "wb");
	if (!file) {
		fprintf(stderr, "Unable to open file %s\n", name);
		return;
	}
	
	//Get HTK header
	header.nSamples = nSamples;
	header.sampPeriod = sampPeriod;
	header.sampSize = sampSize;
	header.parmKind = parmKind;
  #ifdef BIT64
	if (SWAP) {SwapInt(&header.nSamples);} //printf("ns=%i\n", (*header).nSamples);
	if (SWAP) {SwapInt(&header.sampPeriod);} //printf("sp=%i\n", (*header).sampPeriod);
  #else
	if (SWAP) {SwapLong32(&header.nSamples);} //printf("ns=%i\n", header.nSamples);
	if (SWAP) {SwapLong32(&header.sampPeriod);} //printf("sp=%i\n", header.sampPeriod);
	#endif
	if (SWAP) {SwapShort(&header.sampSize);} //printf("ss=%hi\n", header.sampSize);
	if (SWAP) {SwapShort(&header.parmKind);} //printf("pk=%hi\n", header.parmKind);
	
	//Write header
	fwrite(&header, sizeof(htk_header_t), 1, file);
	
	//Write data
	fileLen = sampSize*nSamples; //printf("fileLen=%i\n", fileLen);
	totalparam = fileLen/(sizeof(float));
	buffer = (float *)malloc(fileLen);
        assert(buffer!=NULL);
	for(i = 0; i < totalparam; ++i) {
		if (SWAP) { SwapFloat((data+i)); }
		buffer[i] = *(data+i);
	}
	fwrite(buffer, sizeof(float), totalparam, file);
	
	fclose(file);
	free(buffer);
}

int isDir(const char *target)
{
	#ifdef MACOS //mac (use gnu coreutils tool test, need to call shellscript)
	FILE *fp;
	char b[1];
	char command[strlen(target)+41];
	sprintf(command,"if test -d %s; then echo 1; else echo 0; fi",target);
	fp = popen(command,"r");
	fscanf(fp,"%s",b); /* read output from command */
	fclose(fp);
	return atoi(b);
	#else //linux (need sys/stat.h)
 	struct stat statbuf;
 	if (stat(target, &statbuf)==-1) {printf("stat() error, file \"%s\" not found?\n",target); return -1;}
 	return S_ISDIR(statbuf.st_mode);
	//windows, use FindFirstFile? OR BOOL PathIsDirectory(__in  LPCTSTR pszPath);
	#endif
}

void listdir(char *name, dir_t *directo)
{
	DIR           *d;
	struct dirent *dir;
	int n,numfiles=0;

	chdir(name); // need to change dir to run stat() in isDir()
	d = opendir(name);
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{//printf("--%s %i\n", dir->d_name,isDir(dir->d_name));
			if ((strcmp(dir->d_name, ".")==0) || (strcmp(dir->d_name, "..")==0)) continue;
			if (isDir(dir->d_name)>0) continue;  /* it is a directory */
			numfiles++;
		}
		closedir(d);//printf("closedir:%i\n",closedir(d));
	}

	/* memory allocation */
	(*directo).filenam = (char **)malloc(numfiles * sizeof(char *));
        assert((*directo).filenam!=NULL);

	numfiles=0;
	d = opendir(name);
	if (d)
	{
		while ((dir=readdir(d)) != NULL)
		{
			if ((strcmp(dir->d_name, ".")==0) || (strcmp(dir->d_name, "..")==0)) continue;
			if (isDir(dir->d_name)>0) continue;  /* it is a directory */
			*((*directo).filenam+numfiles) = (char *)malloc((strlen(dir->d_name)+1) * sizeof(char));
                        assert(*((*directo).filenam+numfiles)!=NULL);
			memcpy( *((*directo).filenam+numfiles), dir->d_name, strlen(dir->d_name)+1);
			//printf("__%s\n",*(directo.filenam+numfiles));
			numfiles++;
		}
		closedir(d);//printf("closedir:%i\n",closedir(d));//printf(" HERE ");
	}

	printf("%i files\n", numfiles);
	(*directo).numfiles = numfiles;
}

void UTMatIndex(int N, int idx, int *n, int *m)
{/* convert vector index idx to strickly upper triangular matrix of size N indexes */
	*n = N-2-floor((sqrt(8*( N*(N-1)/2-1-idx)+1)-1)/2);
	*m = idx + (*n+1) - ((2*N-(*n+1))*(*n))/2;
}

int VectUIMatIndex(int N, int n, int m)
{/* convert the coordinate of a strickly upper triangular matrix of size N to vector index */
	return ((2*N-(n+1))*(n))/2 - (n+1) + m;
}

int npairs(int nitem)
{/* number of pairwise comparisons to make (number of combinations)
	int n;
	long double fact1=1, fact2=1; printf("Nitem:%i\n",nitem); 
	for (n=1;n<=nitem;n++) {fact1=fact1*n;} printf("Fact1:%Lf\n",fact1); 
	for (n=1;n<=nitem-2;n++) {fact2=fact2*n;} printf("Fact2:%Lf\n",fact2);
	printf("Npairs:%i\n",(int)(fact1/(2*fact2))); 
	return (int)(fact1/(2*fact2)); */
	return (nitem*(nitem-1)/2);
}

void Align(align_t *alignment, dir_t directo, int id1, int id2, int ce)
{
    int numv1, numv2, n;

    /* memory allocation */
    (*alignment).header1 = (htk_header_t *)malloc(sizeof(htk_header_t));
    assert((*alignment).header1!=NULL);
    (*alignment).header2 = (htk_header_t *)malloc(sizeof(htk_header_t));
    assert((*alignment).header2!=NULL);

    /* read the files */
    (*alignment).mat1 = ReadHTKFile(*(directo.filenam+id1),&(*alignment).numrows,&numv1,(*alignment).header1);
    (*alignment).mat2 = ReadHTKFile(*(directo.filenam+id2),&(*alignment).numcols,&numv2,(*alignment).header2);

    /* controls */
    if (numv1!=numv2) {fprintf(stderr, "Incompatible file formats (different number of parameters)\n");}
    else {(*alignment).numv=numv1;/* Number of rows */}
    #ifdef DEBUG
    printf("#sample1=%i #sample2=%i #parameter=%i\n",(*alignment).numrows,(*alignment).numcols,(*alignment).numv);
    #endif

    /* memory allocation */
    (*alignment).distance = (float *)malloc(sizeof(float));
    assert((*alignment).distance!=NULL);
    (*alignment).cow = (float *)malloc((*alignment).numv * sizeof(float));
    assert((*alignment).cow!=NULL);
    (*alignment).d = (float *)malloc((*alignment).numrows * (*alignment).numcols * sizeof(float));
    assert((*alignment).d!=NULL);
    (*alignment).traceb = (float *)malloc(((*alignment).numrows+1) * ((*alignment).numcols+1) * sizeof(float));
    assert((*alignment).traceb!=NULL);
    (*alignment).indelc = (float *)malloc(sizeof(float));
    assert((*alignment).indelc!=NULL);
    (*alignment).duree = (int *)malloc(sizeof(int));
    assert((*alignment).duree!=NULL);
    
    /* record allocated addresses */
    *(ADDRESSES_HTK_HEADER_T + LASTADD_HTK_HEADER_T++) = (*alignment).header1;
    *(ADDRESSES_HTK_HEADER_T + LASTADD_HTK_HEADER_T++) = (*alignment).header2;
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignment).distance;
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignment).cow;
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignment).d;
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignment).traceb;
    *(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignment).indelc;
    *(ADDRESSES_INT + LASTADD_INT++) = (*alignment).duree;

    /* default vector for parameter factors */
    for (n=0;n<(*alignment).numv;n++) {*((*alignment).cow+n)=1;}

    /* compute pairwise vector distance */
    *(*alignment).indelc = 0 ;
    euclidist((*alignment).d,(*alignment).mat1,(*alignment).mat2,(*alignment).indelc,(*alignment).numrows,(*alignment).numcols,(*alignment).numv,(*alignment).cow);
    /* Call the warpav subroutine */
    if (ce==1) {
      *(*alignment).distance = warpav_ce((*alignment).traceb,(*alignment).d,(*alignment).numrows,(*alignment).numcols,(*alignment).indelc);
    } else {
      *(*alignment).distance = warpav((*alignment).traceb,(*alignment).d,(*alignment).numrows,(*alignment).numcols,(*alignment).indelc);
    }
    /* Call the traceback subroutine to get the average sequence */
    if ((*alignment).weight!=0) {
   	(*alignment).mataverage = averseq((*alignment).duree,(*alignment).traceb,(*alignment).mat1,(*alignment).mat2,(*alignment).weight,(*alignment).numrows,(*alignment).numcols,(*alignment).numv);
		if ( strlen((*alignment).nameave)>0 ) {
			WriteHTKFile((*alignment).nameave, (*alignment).mataverage, (long)*(*alignment).duree, (*(*alignment).header1).sampPeriod, (*(*alignment).header1).sampSize, (*(*alignment).header1).parmKind);
		}
    }

    #ifdef DEBUG
    printf("DISTANCE=%f\n",*((*alignment).distance));
    #else
    printf("%f\n",*((*alignment).distance));
    #endif
}

void free_align_t(align_t *alignment)
{/* free all the fields of a structure align_t but not the pointer to the structure itself */
	free((*alignment).header1); (*alignment).header1 = NULL; /* to avoid Invalid free() later */
	free((*alignment).header2); (*alignment).header2 = NULL;
//printf("\nfree[%p]",(*alignment).mat1);
	//free((*alignment).mat1); (*alignment).mat1 = NULL;
//printf("\nfree[%p]",(*alignment).mat2);
	//free((*alignment).mat2); (*alignment).mat2 = NULL;
	free((*alignment).distance); (*alignment).distance = NULL;
	free((*alignment).cow); (*alignment).cow = NULL;
	free((*alignment).d); (*alignment).d = NULL;
	free((*alignment).traceb); (*alignment).traceb = NULL;
	free((*alignment).indelc); (*alignment).indelc = NULL;
	//if ((*alignment).weight!=0) {
//printf("\nfree[%p]",(*alignment).mataverage);
	    //free((*alignment).mataverage); (*alignment).mataverage = NULL;}
	free((*alignment).duree); (*alignment).duree = NULL;
	free((*alignment).nameave); (*alignment).nameave = NULL;
}

void null_align_t(align_t *alignment)
{/* set to NULL all the ponters of a structure align_t without freeing the memory, to avoid Invalid free() */
	(*alignment).header1 = NULL;
	(*alignment).header2 = NULL;
	//(*alignment).mat1 = NULL;
	//(*alignment).mat2 = NULL;
	(*alignment).distance = NULL;
	(*alignment).cow = NULL;
	(*alignment).d = NULL;
	(*alignment).traceb = NULL;
	(*alignment).indelc = NULL;
	//(*alignment).mataverage = NULL;
	(*alignment).duree = NULL;
	(*alignment).nameave = NULL;
}

void print_align_t(align_t *alignment)
{/* for debugging */
	printf("\n-----");
	printf("\nheader1.nSamples=%li", (long int)(*(*alignment).header1).nSamples);printf(" [%p]",(*alignment).header1);
	printf("\nheader1.sampPeriod=%li", (long int)(*(*alignment).header1).sampPeriod);
 	printf("\nheader1.sampSize=%hi", (*(*alignment).header1).sampSize);
 	printf("\nheader1.parmKind=%hi", (*(*alignment).header1).parmKind);
	printf("\nheader2.nSamples=%li", (long int)(*(*alignment).header2).nSamples);printf(" [%p]",(*alignment).header2);
	printf("\nheader2.sampPeriod=%li", (long int)(*(*alignment).header2).sampPeriod);
 	printf("\nheader2.sampSize=%hi", (*(*alignment).header2).sampSize);
 	printf("\nheader2.parmKind=%hi", (*(*alignment).header2).parmKind);
	printf("\nmat1=%f %f %f ...",*(*alignment).mat1,*((*alignment).mat1+1),*((*alignment).mat1+2));printf(" [%p]",(*alignment).mat1);
	printf("\nmat2=%f %f %f ...",*(*alignment).mat2,*((*alignment).mat2+1),*((*alignment).mat2+2));printf(" [%p]",(*alignment).mat2);
	printf("\ndistance=%f",*(*alignment).distance);printf(" [%p]",(*alignment).distance);
	printf("\ncow=%f %f %f ...",*(*alignment).cow,*((*alignment).cow+1),*((*alignment).cow+2));printf(" [%p]",(*alignment).cow);
	printf("\nd=%f %f %f ...",*(*alignment).d,*((*alignment).d+1),*((*alignment).d+2));printf(" [%p]",(*alignment).d);
	printf("\ntraceb=%f %f %f ...",*(*alignment).traceb,*((*alignment).traceb+1),*((*alignment).traceb+2));printf(" [%p]",(*alignment).traceb);
	printf("\nindelc=%f",*(*alignment).indelc);printf(" [%p]",(*alignment).indelc);
	printf("\nmataverage=%f %f %f ...",*(*alignment).mataverage,*((*alignment).mataverage+1),*((*alignment).mataverage+2));printf(" [%p]",(*alignment).mataverage);
	printf("\nweight=%f",(*alignment).weight);
	printf("\nduree=%i",*(*alignment).duree);printf(" [%p]",(*alignment).duree);
	printf("\nnumrows=%i",(*alignment).numrows);
	printf("\nnumcols=%i",(*alignment).numcols);
	printf("\nnumv=%i",(*alignment).numv);
	printf("\nnameave=%s",(*alignment).nameave);printf(" [%p]",(*alignment).nameave);
	printf("\n-----");
}

int main(int argc, char **argv)
{
    int m, n, i, c, idx, l;
    int npc, npcd, nseq, largest, previous, mode=0, maxal=0, do_ave=0, do_ce=0;
    int *a, *b, *x, *y, *to_be_updated;
    float weight, *pdist;
    dir_t *directo;
    align_t *alignment;
    align_t **alignments;
    FILE *fp;

    /* memory allocation */
    directo = (dir_t *)malloc(sizeof(dir_t));
    assert(directo!=NULL);
    alignment = (align_t *)malloc(sizeof(align_t));
    assert(alignment!=NULL);

    /* default values */
    weight = 0;

    /* arguments check */
    while (1) {
		  static struct option long_options[] = {
		      {"free-end", no_argument, 0, 'f'},
		      {"vector-distance", required_argument, 0, 'd'},
		      {"weight", required_argument, 0, 'w'},
		      {"average", no_argument, 0, 'a'},
		      {"compression-expansion", no_argument, 0, 'c'},
		      {0, 0, 0, 0}
		  };
		  /* getopt_long stores the option index here. */
		  int option_index = 0;
		  c = getopt_long(argc, argv, "fd:w:ac", long_options, &option_index);
		  /* Detect the end of the options. */
		  if (c == -1) break;
		  switch (c) {
			  case 0:
				  /* If this option set a flag, do nothing else now. */
				  if (long_options[option_index].flag != 0)
					  break;
				  printf("option %s", long_options[option_index].name);
				  if (optarg)
					  printf (" with arg %s", optarg);
				  printf("\n");
				  break;
			  case 'f':
				  puts("option -f\n");
				  break;
			  case 'd':
				  printf("option -d with value `%s'\n", optarg);
				  break;
			  case 'w':
				  weight = strtod(optarg,0);
				  break;
			  case 'a':
				  do_ave = 1;
				  break;
			  case 'c':
				  do_ce = 1;
				  break;
			  case '?':/* getopt_long already printed an error message. */
				  break;
			  default:
				  abort();
		  }
    }

    /* check endianality */
    SWAP = CheckArchEndianality(); 
    #ifdef DEBUG
    printf("endianality:%i\n",SWAP);
    printf("optind=%i\n",optind);
    for (n=0;n<=optind;n++) {printf("opt%i=%s\n",n,argv[n]);}
    #endif
    chdir(argv[optind]); // need to change dir to run stat() in isDir()
    if (isDir(argv[optind])) {
    	//compute the average of the files in the directory OR PAIRWISE DISTANCE MATRIX!!!
		mode=1; //to remember for freeing memory
		a = (int *)malloc(sizeof(int));
                assert(a!=NULL);
		b = (int *)malloc(sizeof(int));
                assert(b!=NULL);
		x = (int *)malloc(sizeof(int));
                assert(x!=NULL);
		y = (int *)malloc(sizeof(int));
                assert(y!=NULL);
		/* read directory */
		listdir(argv[optind], directo);
		#ifdef DEBUG
		for (n=0;n<(*directo).numfiles;n++) {printf("--%i %s\n",n,*((*directo).filenam+n));}
		#endif
		npc = npairs((*directo).numfiles);
		#ifdef DEBUG
		printf("%i combinations\n",npc);
		#endif
		pdist = (float *)malloc(npc * sizeof(float));
		assert(pdist!=NULL);
		alignments = (align_t **)malloc(npc * sizeof(align_t));
                assert(alignments!=NULL);
		to_be_updated = (int *)malloc(npc * sizeof(int));
                assert(to_be_updated!=NULL);
		//store memory addresses for mat1, mat2, mataverage, need to calculate the maximum number of alignment possibly required because for each one, a new mataverage is allocated
		for (n=(*directo).numfiles;n>1;n--) {
			maxal = maxal + npairs(n);
		}
		ADDRESSES_HTK_HEADER_T = (htk_header_t **)malloc(2 * maxal * sizeof(htk_header_t *));
		assert(ADDRESSES_HTK_HEADER_T!=NULL);
		ADDRESSES_FLOAT = (float **)malloc(8 * maxal * sizeof(float *));
		assert(ADDRESSES_FLOAT!=NULL);
		ADDRESSES_INT = (int **)malloc(maxal * sizeof(int *));
		assert(ADDRESSES_INT!=NULL);
		ADDRESSES_CHAR = (char **)malloc(maxal * sizeof(char *));
		assert(ADDRESSES_CHAR!=NULL);
		idx = 0;
		fp = fopen("pairwise_distances.csv", "w+");
		for (n=0;n<(*directo).numfiles;n++) {
			fprintf(fp, "%s, ",*((*directo).filenam+n));
			for (l=-1;l<n;l++) {fprintf(fp, ", ");}
			for (m=n+1;m<(*directo).numfiles;m++) {
				#ifdef DEBUG
				printf("idx:%i\n",idx);
				#endif
				//strickly upper triangular matrix may be more compactly represented in an array by storing entry (i,j) in location n(i-1) + j - i(i-1)/2 - i = (2n-i)(i-1)/2 - i + j
				//idx = ((2*((*directo).numfiles)-(n+1))*(n))/2 - (n+1) + m; although it is cheaper to simply increment it
				//UTMatIndex((*directo).numfiles, idx, a, b);
				//printf("%i,%i,%i %i,%i\n",n,m,idx,*a,*b);
				alignments[idx] = (align_t *)malloc(sizeof(align_t));
				assert(alignments[idx]!=NULL);
				(*alignments[idx]).weight = .5;
				(*alignments[idx]).nameave = (char *)malloc(sizeof(char));
				assert((*alignments[idx]).nameave!=NULL);
				*(ADDRESSES_CHAR + LASTADD_CHAR++) = (*alignments[idx]).nameave;
				(*alignments[idx]).nameave[0] = '\0';
				Align(alignments[idx], *directo, n, m, do_ce);
				*(pdist+idx) = *((*alignments[idx]).distance);
				#ifdef DEBUG
				print_align_t(alignments[idx]);
				#endif
				fprintf(fp, "%f, ",*(pdist+idx));
				if (do_ave==0) { //free memory, if no need for average computation
					free(*(ADDRESSES_HTK_HEADER_T + --LASTADD_HTK_HEADER_T));
					free(*(ADDRESSES_HTK_HEADER_T + --LASTADD_HTK_HEADER_T));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_FLOAT + --LASTADD_FLOAT));
					free(*(ADDRESSES_INT + --LASTADD_INT));
					free(*(ADDRESSES_CHAR + --LASTADD_CHAR));
				}
				idx++;
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		#ifdef DEBUG
		printf("%i alignments\n",idx);
		for (n=0;n<npc;n++) {printf(" %f",*(pdist+n));}
		#endif
		printf("\n");
		if (do_ave==1) {
			nseq = (*directo).numfiles;
			npcd = npc; //need npc for freeing memory so we duplicate it here!
			while (nseq>1) {/* recursively compute average sequences, until only 2 sequences remain
				a is always inferior to b so the largest average gets rank a in alignments and pdist */
				largest = maxcost(pdist,0,npcd-1);
				UTMatIndex(nseq,largest,a,b);
				npcd = npairs(--nseq);
				#ifdef DEBUG
				printf("\n\nmaxi=%i (%i,%i)",largest,*a,*b);
				printf("\nnseq=%i",nseq);
				#endif
				/* update all the mat1 and mat2 */
				for (n=0;n<npcd;n++) {
					/* the average sequence from the largest distance alignment replaces mat1 in this alignment structure
					thus, it takes the index of the sequence mat1 and we need to update the alignments involving it */
					UTMatIndex(nseq,n,x,y); //printf("\n-> n:%i, *x:%i, *y:%i",n,*x,*y);
					if (*x==*a) { /* need to perform a new alignment */
						to_be_updated[n] = 1;
						/* use the average sequence as a new sequence */
						(*alignments[n]).numrows = *(*alignments[largest]).duree; //printf("\nduree=%i",(*alignments[n]).numcols);
						(*alignments[n]).mat1 = (*alignments[largest]).mataverage; 
						/* find a previous alignment with the length and parameters of the second sequence */
						*y = (*y<=*b) ? *y : *y+1 ; //need to increment because of index change
						previous = VectUIMatIndex(nseq+1, *a, *y);
						if (n!=previous) {
							#ifdef DEBUG
							printf("\n(*x==*a) align average of [%i,%i](al.idx:%i) with mat2 from alignment idx:%i(%i,%i)",*a,*b,largest,previous,*a,*y);
							#endif
							(*alignments[n]).numcols = (*alignments[previous]).numcols;
							(*alignments[n]).mat2 = (*alignments[previous]).mat2;
						}	
						#ifdef DEBUG				
						print_align_t(alignments[n]);
						#endif
					} else if (*y==*a) { /* need to perform a new alignment */
						to_be_updated[n] = 1;
						/* use the average sequence as a new sequence */
						(*alignments[n]).numcols = *(*alignments[largest]).duree;
						(*alignments[n]).mat2 = (*alignments[largest]).mataverage;
						/* find a previous alignment with the length and parameters of the second sequence */
						*x = (*x<=*b) ? *x : *x+1 ; //need to increment because of index change
						previous = VectUIMatIndex(nseq+1, *x, *a);
						if (n!=previous) {
							#ifdef DEBUG
							printf("\n(*y==*a) align average of [%i,%i](al.idx:%i) with mat1 from alignment idx:%i(%i,%i)",*a,*b,largest,previous,*x,*a);
							#endif
							(*alignments[n]).numrows = (*alignments[previous]).numrows;
							(*alignments[n]).mat1 = (*alignments[previous]).mat1;
						}
					} else { /* copy the proper alignments with change of index */
						to_be_updated[n] = 0;
						*x = (*x<*b) ? *x : *x+1;
						*y = (*y<*b) ? *y : *y+1;
						previous = VectUIMatIndex(nseq+1, *x, *y);
						if (n!=previous) {
							#ifdef DEBUG
							printf("\n(*x!=*a && *y!=*a) (*a=%i,*b=%i) copy alignment idx:%i(%i,%i) to idx:%i",*a,*b,previous,*x,*y,n);
							#endif
							(*alignments[n]) = (*alignments[previous]);
							*(pdist+n) = *((*alignments[n]).distance);
						}
					}
				}
				/* update the mataverage */
				for (n=0;n<npcd;n++) {
					if (to_be_updated[n]==1) {
						//printf("\nALIGNMENT UPDATE %i",n);
						(*alignments[n]).d = (float *)malloc(((*alignments[n]).numrows) * ((*alignments[n]).numcols) * sizeof(float));
						assert((*alignments[n]).d!=NULL);
						*(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignments[n]).d;
						(*alignments[n]).traceb = (float *)malloc(((*alignments[n]).numrows+1) * ((*alignments[n]).numcols+1) * sizeof(float));
						assert((*alignments[n]).traceb!=NULL);
						*(ADDRESSES_FLOAT + LASTADD_FLOAT++) = (*alignments[n]).traceb;
						*(*alignments[n]).indelc = 0;
						euclidist((*alignments[n]).d,(*alignments[n]).mat1,(*alignments[n]).mat2,(*alignments[n]).indelc,(*alignments[n]).numrows,(*alignments[n]).numcols,(*alignments[n]).numv,(*alignments[n]).cow);
						if (do_ce==1) {
						  *(*alignments[n]).distance = warpav_ce((*alignments[n]).traceb,(*alignments[n]).d,(*alignments[n]).numrows,(*alignments[n]).numcols,(*alignments[n]).indelc);
						} else {
						  *(*alignments[n]).distance = warpav((*alignments[n]).traceb,(*alignments[n]).d,(*alignments[n]).numrows,(*alignments[n]).numcols,(*alignments[n]).indelc);
						}
						(*alignments[n]).mataverage = averseq((*alignments[n]).duree,(*alignments[n]).traceb,(*alignments[n]).mat1,(*alignments[n]).mat2,(*alignments[n]).weight,(*alignments[n]).numrows,(*alignments[n]).numcols,(*alignments[n]).numv); //printf("\naverseq duree=%i",*(*alignments[n]).duree);
						*(pdist+n) = *((*alignments[n]).distance);
						to_be_updated[n]=0;
					}
				}
				#ifdef DEBUG
				printf("\n");
				for (n=0;n<npcd;n++) {printf(" %f",*(pdist+n));}
				#endif
			}
			WriteHTKFile("average", (*alignments[0]).mataverage, (long)*(*alignments[0]).duree, (*(*alignments[0]).header1).sampPeriod, (*(*alignments[0]).header1).sampSize, (*(*alignments[0]).header1).parmKind);
		}
    } else {
		(*alignment).weight = weight;
		(*directo).numfiles = 2;
		(*directo).filenam = (char **)malloc(2 * sizeof(char *));
		assert((*directo).filenam!=NULL);
		for (n=0;n<(*directo).numfiles;n++) { 
			*((*directo).filenam+n) = (char *)malloc((strlen(argv[optind+n])+1) * sizeof(char));
			assert(*((*directo).filenam+n)!=NULL);
		}
		memcpy( *((*directo).filenam), argv[optind], strlen(argv[optind])+1 );
		memcpy( *((*directo).filenam+1), argv[optind+1], strlen(argv[optind+1])+1 );
		ADDRESSES_HTK_HEADER_T = (htk_header_t **)malloc(2 * 2 * sizeof(htk_header_t *));
		assert(ADDRESSES_HTK_HEADER_T!=NULL);
		ADDRESSES_FLOAT = (float **)malloc(8 * 2 * sizeof(float *));
		assert(ADDRESSES_FLOAT!=NULL);
		ADDRESSES_INT = (int **)malloc(2 * sizeof(int *));
		assert(ADDRESSES_INT!=NULL);
		ADDRESSES_CHAR = (char **)malloc(2 * sizeof(char *));
		assert(ADDRESSES_CHAR!=NULL);
		if (argv[optind+2]!=NULL) {//compute and write an average sequence
			(*alignment).nameave = (char *)malloc((strlen(argv[optind+2])+1) * sizeof(char));
		        assert((*alignment).nameave!=NULL);
			*(ADDRESSES_CHAR + LASTADD_CHAR++) = (*alignment).nameave;
			memcpy( (*alignment).nameave, argv[optind+2], strlen(argv[optind+2])+1 );
			Align(alignment, *directo, 0, 1, do_ce);
		} else {//only compute the distance
			//for (n=0;n<(*directo).numfiles;n++) {printf("--%s\n",*((*directo).filenam+n));}
			(*alignment).weight = 0;
			(*alignment).nameave = (char *)malloc(sizeof(char));
                        assert((*alignment).nameave!=NULL);
			*(ADDRESSES_CHAR + LASTADD_CHAR++) = (*alignment).nameave;
			(*alignment).nameave[0] = '\0';
			Align(alignment, *directo, 0, 1, do_ce);
		}
    }

	if (mode==1) {
		free(a);
		free(b);
		free(x);
		free(y);
		free(pdist);
		for (n=0;n<npc;n++) {
			free(alignments[n]);
		}
		free(alignments);
		free(to_be_updated);
	}
	free(alignment);
	for (n=0;n<(*directo).numfiles;n++) {
		free(*((*directo).filenam+n));
	}
	free((*directo).filenam);
	free(directo);
	for (n=0;n<LASTADD_HTK_HEADER_T;n++) {
		free(*(ADDRESSES_HTK_HEADER_T+n));
	}
	free(ADDRESSES_HTK_HEADER_T);
	for (n=0;n<LASTADD_FLOAT;n++) {
		free(*(ADDRESSES_FLOAT+n));
	}
	free(ADDRESSES_FLOAT);
	for (n=0;n<LASTADD_INT;n++) {
		free(*(ADDRESSES_INT+n));
	}
	free(ADDRESSES_INT);
	for (n=0;n<LASTADD_CHAR;n++) {
		free(*(ADDRESSES_CHAR+n));
	}
	free(ADDRESSES_CHAR);

    return 1 ;
}


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "KRnorm.h"


int main(int argc, char *argv[])
{
    int resolution;
    FILE* gsize = NULL;
    FILE* digest = NULL;
    FILE* contact = NULL;

    FILE* norm1=NULL;
    FILE* norm2=NULL;

    int inter;
    long long chr1_size = 0;
    long long chr2_size = 0;
    long long size;
    char chr1[64];
    char chr2[64];
    char chr[64];
    int frag1;
    int frag2;

    int* chr1_starts = NULL;
    int* chr1_ends   = NULL;
    int* chr2_starts = NULL;
    int* chr2_ends   = NULL;
    int  chr1_base   = 0;
    int  chr2_base   = 0;

    int  chr1_frag_num = 0;
    int  chr2_frag_num = 0;

    int start, start1, start2;
    int end, end1, end2;
    int bin1_start, bin1_end, bin2_start, bin2_end;
    int id, i, j;

    float* matrix = NULL;
    float* norm   = NULL;


    if (argc != 5 && argc != 6 && argc != 7)
    {
        fprintf(stderr, "Usage: %s <resolution> <genome size> <digest file> <contact pairs> [<norm out file>]\n\n", argv[0]);
        return 1;
    }
    sscanf(argv[1], "%d", &resolution);

    if ((gsize=fopen(argv[2], "r")) == NULL)
    {
        fprintf(stderr, "Cannot open genome size file: %s\n\n", argv[2]);
        return 1;
    }
    if ((digest=fopen(argv[3], "r")) == NULL)
    {
        fprintf(stderr, "Cannot open genome digest file: %s\n\n", argv[3]);
        return 1;
    }
    if ((contact=fopen(argv[4], "r")) == NULL)
    {
        fprintf(stderr, "Cannot open contact pairs file: %s\n\n", argv[4]);
        return 1;
    }
    fscanf(contact, "%s %d %s %d", chr1, &frag1, chr2, &frag2);
    rewind(contact);
    inter = strcmp(chr1, chr2);
    if (argc > 5 && (norm1=fopen(argv[5], "w")) == NULL)
    {
        fprintf(stderr, "Cannot create normalizing file %s for %s\n\n", argv[5], chr1);
        return 1;
    }
    if (inter && argc == 7 && (norm2=fopen(argv[6], "w")) == NULL)
    {
        fprintf(stderr, "Cannot create normalizing file %s for %s\n\n", argv[6], chr2);
        return 1;
    }


    while (fscanf(gsize, "%s %lld", chr, &size) != EOF)
    {
        if (strcmp(chr, chr1) == 0) chr1_size = size;
        if (strcmp(chr, chr2) == 0) chr2_size = size;
        if (chr1_size && chr2_size) break;
    }
    fclose(gsize);

    size = chr1_size / resolution;
    if (chr1_size % resolution) size += 1;
    chr1_size = size;
    size = chr2_size / resolution;
    if (chr2_size % resolution) size += 1;
    chr2_size = size;
    if (inter)
        size = chr1_size + chr2_size;
    else
        size = chr1_size;


    while (fscanf(digest, "%s %d %d %d", chr, &start, &end, &id) != EOF)
    {
        if (strcmp(chr, chr1) == 0)
        {
            if (chr1_base == 0)
                chr1_base = id;
            chr1_frag_num++;
        }
        else
            if (strcmp(chr, chr2) == 0)
            {
                if (chr2_base == 0)
                    chr2_base = id;
                chr2_frag_num++;
            }
            else
                if ((chr2_base != 0) || (chr1_base != 0 && inter == 0))
                    break;
    }
    rewind(digest);

    if (chr1_frag_num)
    {
        chr1_starts = (int *)malloc(chr1_frag_num * sizeof(int));
        chr1_ends   = (int *)malloc(chr1_frag_num * sizeof(int));
    }
    if (chr2_frag_num)
    {
        chr2_starts = (int *)malloc(chr2_frag_num * sizeof(int));
        chr2_ends   = (int *)malloc(chr2_frag_num * sizeof(int));
    }
    else
    {
        chr2_starts = chr1_starts;
        chr2_ends   = chr1_ends;
        chr2_base   = chr1_base;
    }


    while (fscanf(digest, "%s %d %d %d", chr, &start, &end, &id) != EOF)
    {
        if (strcmp(chr, chr1) == 0)
        {
            chr1_starts[id - chr1_base] = start;
            chr1_ends[id - chr1_base] = end;
            chr1_frag_num--;
        }
        else
            if (strcmp(chr, chr2) == 0)
            {
                chr2_starts[id - chr2_base] = start;
                chr2_ends[id - chr2_base] = end;
                chr2_frag_num--;
            }
            else
                if (chr1_frag_num == 0 && chr2_frag_num == 0)
                    break;
    }
    fclose(digest);


    matrix = (float *)calloc(size * size, sizeof(float));


    while (fscanf(contact, "%s %d %s %d", chr, &frag1, chr, &frag2) != EOF)
    {
        start1 = chr1_starts[frag1 - chr1_base];
        end1   = chr1_ends[frag1 - chr1_base] - 1;
        start2 = chr2_starts[frag2 - chr2_base];
        end2   = chr2_ends[frag2 - chr2_base] - 1;
        bin1_start = start1 / resolution;
        bin1_end   = end1 / resolution;
        bin2_start = start2 / resolution;
        bin2_end   = end2 / resolution;
        if (inter)
        {
            bin2_start += chr1_size;
            bin2_end   += chr1_size;
        }

        matrix[size * bin1_start + bin2_start]++;
        matrix[size * bin2_start + bin1_start]++;
        if (bin2_end != bin2_start)
        {
            matrix[size * bin1_start + bin2_end]++;
            matrix[size * bin2_end + bin1_start]++;
        }
        if (bin1_end != bin1_start)
        {
            matrix[size * bin1_end + bin2_start]++;
            matrix[size * bin2_start + bin1_end]++;
            if (bin2_end != bin2_start)
            {
                matrix[size * bin1_end + bin2_end]++;
                matrix[size * bin2_end + bin1_end]++;
            }
        }
    }
    fclose(contact);

    if (chr1_starts != NULL) free(chr1_starts);
    if (chr1_ends   != NULL) free(chr1_ends);
    if (inter)
    {
        if (chr2_starts != NULL) free(chr2_starts);
        if (chr2_ends   != NULL) free(chr2_ends);
    }

    //Output
    for (i=0; i<chr1_size; i++)
        for (j=(inter? chr1_size : (i+1)); j<size; j++)
            if (matrix[size * i + j] > 0)
            {
                start = i * resolution;
                end   = (inter? j-chr1_size : j) * resolution;
                printf("%d\t%d\t%.0f\n", start, end, matrix[size * i + j]);
            }

    if (argc > 5)
    {
        for (i=0; i<size; i++)
            matrix[size * i + i] = 1.0;

        norm = (float *)calloc(size, sizeof(float));

        KRnorm(matrix, size, norm);
    }

    if (matrix != NULL) free(matrix);


    if (argc > 5)
    {
        for (i=0; i<chr1_size; i++)
            fprintf(norm1, "%f\n", 1.0 / norm[i]);
        fclose(norm1);
        if (inter)
        {
            for (i=chr1_size; i<size; i++)
                fprintf(norm2, "%f\n", 1.0 / norm[i]);
            fclose(norm2);
        }
    }

    if (norm != NULL) free(norm);

    return 0;
}




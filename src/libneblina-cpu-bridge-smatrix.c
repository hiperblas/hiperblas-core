#include "libneblina.h"
#include <stdio.h>
#include <stdlib.h>

slist *slist_add(slist *list, int col, double re, double im)
{
    slist *nlist = (slist *)malloc(sizeof(slist));
    nlist->col = col;
    nlist->re = re;
    nlist->im = im;
    nlist->next = list;
    return nlist;
}

void slist_clear(slist *list)
{
    if (list == NULL)
        return;
    do
    {
        slist *tmp = list->next;
        free(list);
        list = tmp;
    } while (list != NULL);
}

void smatrix_pack(smatrix_t *smatrix)
{
    // essa parte seta os valores que foram lidos do arquivo em smat para o formato de
    // array continuo
    smatrix->idx_col = (int *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(int));
    smatrix->m = (double *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(double));

    memset(smatrix->idx_col, -1, smatrix->nrow * smatrix->maxcols * sizeof(int));
    memset(smatrix->m, 0, smatrix->nrow * smatrix->maxcols * sizeof(double));

    slist *tmp = NULL;
    for (int i = 0; i < smatrix->nrow; i++)
    {
        tmp = smatrix->smat[i];

        while (tmp != NULL)
        {
            int idx = i * smatrix->maxcols + smatrix->icount[i];
            smatrix->idx_col[idx] = tmp->col;
            int midx = i * smatrix->maxcols + smatrix->icount[i];
            smatrix->m[midx] = tmp->re;
            smatrix->icount[i]++;
            tmp = tmp->next;
        }
        slist_clear(smatrix->smat[i]);
        smatrix->smat[i] = NULL;
    }
    smatrix->isPacked = 1;
}

void smatrix_pack_complex(smatrix_t *smatrix)
{
    smatrix->idx_col = (int *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(int));
    smatrix->m = (double *)malloc(smatrix->nrow * smatrix->maxcols * COMPLEX_SIZE);

    memset(smatrix->idx_col, -1, smatrix->nrow * smatrix->maxcols * sizeof(int));
    memset(smatrix->m, 0, smatrix->nrow * smatrix->maxcols * COMPLEX_SIZE);

    slist *tmp = NULL;
    for (int i = 0; i < smatrix->nrow; i++)
    {
        tmp = smatrix->smat[i];
        while (tmp != NULL)
        {
            int idx = i * smatrix->maxcols + smatrix->icount[i];
            smatrix->idx_col[idx] = tmp->col;
            int midx = 2 * (i * smatrix->maxcols + smatrix->icount[i]);
            smatrix->m[midx] = tmp->re;
            smatrix->m[midx + 1] = tmp->im;
            smatrix->icount[i]++;
            tmp = tmp->next;
        }
        slist_clear(smatrix->smat[i]);
        smatrix->smat[i] = NULL;
    }
    smatrix->isPacked = 1;
}

void smatrix_set_real_value(smatrix_t *smatrix, int i, int j, double r)
{
    if (smatrix->smat == NULL)
    {
        smatrix->smat = (slist **)malloc(smatrix->nrow * sizeof(slist *));
        for (int x = 0; x < smatrix->nrow; x++)
        {
            smatrix->smat[x] = NULL;
        }
    }

    if (i < 0 || i >= smatrix->nrow)
    {
        printf("invalid row index on loading sparse matrix %d\n", i);
        exit(-1);
    }
    if (j < 0 || j >= smatrix->ncol)
    {
        printf("invalid col index on loading sparse matrix\n");
        exit(-1);
    }

    smatrix->smat[i] = slist_add(smatrix->smat[i], j, r, 0.0);
    smatrix->rcount[i]++;
    if (smatrix->rcount[i] > smatrix->maxcols)
    {
        smatrix->maxcols = smatrix->rcount[i];
    }
}

void smatrix_set_complex_value(smatrix_t *smatrix, int i, int j, double r, double im)
{
    if (smatrix->smat == NULL)
    {
        smatrix->smat = (slist **)malloc(smatrix->nrow * sizeof(slist *));
        for (int x = 0; x < smatrix->nrow; x++)
            smatrix->smat[x] = NULL;
    }

    if (i < 0 || i >= smatrix->nrow)
    {
        printf("invalid row index on loading sparse matrix %d\n", i);
        exit(-1);
    }
    if (j < 0 || j >= smatrix->ncol)
    {
        printf("invalid col index on loading sparse matrix\n");
        exit(-1);
    }

    smatrix->smat[i] = slist_add(smatrix->smat[i], j, r, im);
    smatrix->rcount[i]++;
    if (smatrix->rcount[i] > smatrix->maxcols)
    {
        smatrix->maxcols = smatrix->rcount[i];
    }
}

smatrix_t *smatrix_new(int nrow, int ncol, data_type type)
{
    smatrix_t *smatrix = (smatrix_t *)malloc(sizeof(smatrix_t));
    smatrix->ncol = ncol;
    smatrix->nrow = nrow;
    smatrix->type = type;

    smatrix->rcount = (int *)malloc(smatrix->nrow * sizeof(int));
    smatrix->icount = (int *)malloc(smatrix->nrow * sizeof(int));

    memset(smatrix->rcount, 0, smatrix->nrow * sizeof(int));
    memset(smatrix->icount, 0, smatrix->nrow * sizeof(int));

    smatrix->location = LOCHOS;
    smatrix->smat = NULL;
    smatrix->extra = NULL;
    smatrix->idxColMem = NULL;
    smatrix->m = NULL;
    smatrix->isPacked = 0;
    smatrix->maxcols = 0;

    return smatrix;
}

void smatrix_t_clear(smatrix_t *smatrix)
{
    if (smatrix != NULL)
    {
        free(smatrix->m);
        free(smatrix->idx_col);
        free(smatrix);
    }
}

void smatrix_load_double(smatrix_t *smatrix, FILE *f)
{
    double e = 0.0;
    int *rcount = (int *)malloc(smatrix->nrow * sizeof(int));
    int *icount = (int *)malloc(smatrix->nrow * sizeof(int));
    slist **smat = (slist **)malloc(smatrix->nrow * sizeof(slist *));
    int i = 0, j = 0;
    ;
    int maxcols = 0;
    memset(rcount, 0, smatrix->nrow * sizeof(int));
    memset(icount, 0, smatrix->nrow * sizeof(int));

    for (i = 0; i < smatrix->nrow; i++)
        smat[i] = NULL;

    while (fscanf(f, "%d %d %lf", &i, &j, &e) != EOF)
    {
        i--;
        j--;
        if (i < 0 || i >= smatrix->nrow)
        {
            printf("invalid row index on loading sparse matrix\n");
            exit(-1);
        }
        if (j < 0 || j >= smatrix->ncol)
        {
            printf("invalid col index on loading sparse matrix\n");
            exit(-1);
        }

        smat[i] = slist_add(smat[i], j, e, 0.0);
        rcount[i]++;
        if (rcount[i] > maxcols)
        {
            maxcols = rcount[i];
        }
    }

    smatrix->maxcols = maxcols;
    smatrix->idx_col = (int *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(int));
    smatrix->m = (double *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(double));

    memset(smatrix->idx_col, -1, smatrix->nrow * smatrix->maxcols * sizeof(int));

    slist *tmp = NULL;
    for (i = 0; i < smatrix->nrow; i++)
    {
        tmp = smat[i];
        while (tmp != NULL)
        {
            smatrix->idx_col[i * maxcols + icount[i]] = tmp->col;
            smatrix->m[i * maxcols + icount[i]] = tmp->re;
            icount[i]++;
            tmp = tmp->next;
        }
        slist_clear(smat[i]);
    }
    free(rcount);
    free(icount);
    free(smat);
    fclose(f);
}

void smatrix_load_complex(smatrix_t *smatrix, FILE *f)
{
    double re, im;
    int *rcount = (int *)malloc(smatrix->nrow * sizeof(int));
    int *icount = (int *)malloc(smatrix->nrow * sizeof(int));
    slist **smat = (slist **)malloc(smatrix->nrow * sizeof(slist *));
    int i = 0, j = 0;
    int maxcols = 0;
    memset(rcount, 0, smatrix->nrow * sizeof(int));
    memset(icount, 0, smatrix->nrow * sizeof(int));

    for (i = 0; i < smatrix->nrow; i++)
    {
        smat[i] = NULL;
    }

    while (fscanf(f, "%d %d %lf %lf", &i, &j, &re, &im) != EOF)
    {
        i--;
        j--;
        if (i < 0 || i >= smatrix->nrow)
        {
            printf("invalid row on loading sparse matrix\n");
            exit(-1);
        }
        if (j < 0 || j >= smatrix->ncol)
        {
            printf("invalid col on loading sparse matrix\n");
            exit(-1);
        }

        smat[i] = slist_add(smat[i], j, re, im);
        rcount[i]++;
        if (rcount[i] > maxcols)
        {
            maxcols = rcount[i];
        }
    }

    smatrix->maxcols = maxcols;
    smatrix->idx_col = (int *)malloc(smatrix->nrow * smatrix->maxcols * sizeof(int));
    smatrix->m = (double *)malloc(2 * smatrix->nrow * smatrix->maxcols * sizeof(double));

    memset(smatrix->idx_col, -1, smatrix->nrow * smatrix->maxcols * sizeof(int));

    slist *tmp = NULL;
    for (i = 0; i < smatrix->nrow; i++)
    {
        tmp = smat[i];
        while (tmp != NULL)
        {
            smatrix->idx_col[i * maxcols + icount[i]] = tmp->col;
            smatrix->m[2 * (i * maxcols + icount[i])] = tmp->re;
            smatrix->m[2 * (i * maxcols + icount[i]) + 1] = tmp->im;
            icount[i]++;
            tmp = tmp->next;
        }
        slist_clear(smat[i]);
    }
    free(rcount);
    free(icount);
    free(smat);
    fclose(f);
}

void smatreqhost(smatrix_t *smatrix)
{
    if (smatrix->location != LOCHOS)
    {
        smatrix->location = LOCHOS;
        smatrix->idx_col = smatrix->idxColMem;
        smatrix->m = smatrix->extra;
        smatrix->extra = NULL;
        smatrix->idxColMem = NULL;
    }
}

void smatreqdev(smatrix_t *smatrix)
{
    if (smatrix->location != LOCDEV)
    {
        smatrix->location = LOCDEV;
        smatrix->idxColMem = smatrix->idx_col;
        smatrix->extra = smatrix->m;
    }
}

void smatrix_delete(smatrix_t *smatrix)
{
    if (smatrix->smat != NULL)
    {
        slist *tmp;
        for (int i = 0; i < smatrix->nrow; i++)
        {
            tmp = smatrix->smat[i];
            slist_clear(smatrix->smat[i]);
            smatrix->smat[i] = NULL;
        }
        free(smatrix->smat);
    }

    if (smatrix->idx_col != NULL)
    {
        free(smatrix->idx_col);
    }

    if (smatrix->rcount != NULL)
    {
        free(smatrix->rcount);
    }

    if (smatrix->icount != NULL)
    {
        free(smatrix->icount);
    }

    if (smatrix->m != NULL)
    {
        free(smatrix->m);
    }
    free(smatrix);
}

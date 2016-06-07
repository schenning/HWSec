/*
 * Copyright (C) Telecom ParisTech
 * 
 * This file must be used under the terms of the CeCILL. This source
 * file is licensed as described in the file COPYING, which you should
 * have received as part of this distribution. The terms are also
 * available at:
 * http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#include <utils.h>
#include <tr_pcc.h>

extern float *tr_new_trace_1 (int l);
extern void tr_free_trace_1 (int l, float *t);
extern void tr_init_trace_1 (int l, float *t, float val);
extern void tr_copy_1 (int l, float *dest, float *src);
extern void tr_acc_1 (int l, float *dest, float *src);
extern void tr_add_1 (int l, float *dest, float *src1, float *src2);
extern void tr_sub_1 (int l, float *dest, float *src1, float *src2);
extern void tr_scalar_mul_1 (int l, float *dest, float *src, float val);
extern void tr_scalar_div_1 (int l, float *dest, float *src, float val);
extern void tr_mul_1 (int l, float *dest, float *src1, float *src2);
extern void tr_div_1 (int l, float *dest, float *src1, float *src2);
extern void tr_sqr_1 (int l, float *dest, float *src);
extern void tr_sqrt_1 (int l, float *dest, float *src);
extern void tr_abs_1 (int l, float *dest, float *src);
extern float tr_min_1 (int l, float *t, int *idx);
extern float tr_max_1 (int l, float *t, int *idx);
extern void tr_print_1 (int l, float *t);
extern void tr_fprint_1 (int l, FILE * fp, float *t);

tr_pcc_context
tr_pcc_init (int l, int ny)
{
  int i;
  tr_pcc_context ctx;

  if (l < 1)
    {
      ERROR (NULL, -1, "Invalid vector length: %d", l);
    }
  if (ny < 1)
    {
      ERROR (NULL, -1, "Invalid number of Y random variables: %d", ny);
    }
  ctx = XCALLOC (1, sizeof (struct tr_pcc_context_s));
  ctx->ny = ny;
  ctx->nr = 0;
  ctx->l = l;
  ctx->rx = tr_new_trace_1 (l);
  ctx->x = tr_new_trace_1 (l);
  tr_init_trace_1 (l, ctx->x, 0.0);
  ctx->x2 = tr_new_trace_1 (l);
  tr_init_trace_1 (l, ctx->x2, 0.0);
  ctx->y = tr_new_trace_1 (ny);
  tr_init_trace_1 (ny, ctx->y, 0.0);
  ctx->y2 = tr_new_trace_1 (ny);
  tr_init_trace_1 (ny, ctx->y2, 0.0);
  ctx->xy = XCALLOC (ny, sizeof (float *));
  ctx->pcc = XCALLOC (ny, sizeof (float *));
  for (i = 0; i < ny; i++)
    {
      ctx->xy[i] = tr_new_trace_1 (l);
      tr_init_trace_1 (l, ctx->xy[i], 0.0);
      ctx->pcc[i] = tr_new_trace_1 (l);
      tr_init_trace_1 (l, ctx->pcc[i], 0.0);
    }
  ctx->state = 0;
  ctx->cnt = ny;
  ctx->flags = XCALLOC (ny, sizeof (char));
  for (i = 0; i < ny; i++)
    {
      ctx->flags[i] = 0;
    }
  ctx->tmp1 = tr_new_trace_1 (l);
  ctx->tmp2 = tr_new_trace_1 (l);
  ctx->tmp3 = tr_new_trace_1 (l);
  return ctx;
}

void
tr_pcc_insert_x (tr_pcc_context ctx, float *x)
{
  if (ctx->cnt != ctx->ny)
    {
      ERROR (, -1, "missing %d Y realizations", ctx->ny - ctx->cnt);
    }
  ctx->cnt = 0;
  ctx->state = 1 - ctx->state;
  tr_copy_1 (ctx->l, ctx->rx, x);
  tr_acc_1 (ctx->l, ctx->x, x);
  tr_sqr_1 (ctx->l, ctx->tmp1, x);
  tr_acc_1 (ctx->l, ctx->x2, ctx->tmp1);
  ctx->nr += 1;
}

void
tr_pcc_insert_y (tr_pcc_context ctx, int ny, float y)
{
  if (ny < 0 || ny >= ctx->ny)
    {
      ERROR (, -1, "Invalid Y index: %d", ny);
    }
  if (ctx->flags[ny] == ctx->state)
    {
      ERROR (, -1, "Y realization #%d inserted twice", ny);
    }
  ctx->y[ny] += y;
  ctx->y2[ny] += y * y;
  tr_scalar_mul_1 (ctx->l, ctx->tmp1, ctx->rx, y);
  tr_acc_1 (ctx->l, ctx->xy[ny], ctx->tmp1);
  ctx->cnt += 1;
  ctx->flags[ny] = ctx->state;
}

void
tr_pcc_consolidate (tr_pcc_context ctx)
{
  float n;
  int i;

  if (ctx->cnt != ctx->ny)
    {
      ERROR (, -1, "missing %d Y realizations", ctx->ny - ctx->cnt);
    }
  if (ctx->nr < 2)
    {
      ERROR (, -1, "not enough realizations (%d, min 2)", ctx->nr);
    }
  n = (float) (ctx->nr);
  tr_scalar_mul_1 (ctx->l, ctx->tmp1, ctx->x2, n);	/* TMP1 = N.X2 */
  tr_sqr_1 (ctx->l, ctx->tmp2, ctx->x);	/* TMP2 = X^2 */
  tr_sub_1 (ctx->l, ctx->tmp1, ctx->tmp1, ctx->tmp2);	/* TMP1 = N.X2-X^2 */
  tr_sqrt_1 (ctx->l, ctx->tmp1, ctx->tmp1);	/* TMP1 = SQRT(N.X2-X^2) */
  for (i = 0; i < ctx->ny; i++)
    {
      tr_scalar_mul_1 (ctx->l, ctx->tmp2, ctx->xy[i], n);	/* TMP2 = N.XY */
      tr_scalar_mul_1 (ctx->l, ctx->tmp3, ctx->x, ctx->y[i]);	/* TMP3 = X.Y */
      tr_sub_1 (ctx->l, ctx->tmp2, ctx->tmp2, ctx->tmp3);	/* TMP2 = N.XY-X.Y */
      tr_div_1 (ctx->l, ctx->tmp2, ctx->tmp2, ctx->tmp1);
      /* TMP2 = (N.XY-X.Y)/(SQRT(N.X2-X^2)) */
      tr_scalar_div_1 (ctx->l, ctx->pcc[i], ctx->tmp2,
		       sqrt (n * ctx->y2[i] - ctx->y[i] * ctx->y[i]));
      /* PCC = (N.XY-X.Y)/(SQRT(N.X2-X^2))/SQRT(N.Y2-Y^2) */
    }
}

float *
tr_pcc_get_pcc (tr_pcc_context ctx, int ny)
{
  if (ny < 0 || ny >= ctx->ny)
    {
      ERROR (NULL, -1, "Invalid Y index: %d", ny);
    }
  return ctx->pcc[ny];
}

void
tr_pcc_free (tr_pcc_context ctx)
{
  int i;

  tr_free_trace_1 (ctx->l, ctx->rx);
  tr_free_trace_1 (ctx->l, ctx->x);
  tr_free_trace_1 (ctx->l, ctx->x2);
  tr_free_trace_1 (ctx->ny, ctx->y);
  tr_free_trace_1 (ctx->ny, ctx->y2);
  for (i = 0; i < ctx->ny; i++)
    {
      tr_free_trace_1 (ctx->l, ctx->xy[i]);
      tr_free_trace_1 (ctx->l, ctx->pcc[i]);
    }
  free (ctx->flags);
  tr_free_trace_1 (ctx->l, ctx->tmp1);
  tr_free_trace_1 (ctx->l, ctx->tmp2);
  tr_free_trace_1 (ctx->l, ctx->tmp3);
  free (ctx);
}

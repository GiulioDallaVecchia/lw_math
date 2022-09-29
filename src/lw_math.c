/******************************************************************************
 * Filename              :   lw_math.c
 * Author                :   Giulio Dalla Vecchia
 * Origin Date           :   5 may 2022
 *
 * Copyright (c) 2022 Giulio Dalla Vecchia. All rights reserved.
 *
 ******************************************************************************/

/** @file lw_math.c
 *  @brief This module handles the basic math operations in fixed-point format
 */

/*****************************************************************************
 * Includes
 ******************************************************************************/
#include "lw_math.h"

/*****************************************************************************
 * Module Preprocessor Constants
 ******************************************************************************/

#define SIN_COS_TABLE {\
0x0000,0x00C9,0x0192,0x025B,0x0324,0x03ED,0x04B6,0x057F,\
0x0648,0x0711,0x07D9,0x08A2,0x096A,0x0A33,0x0AFB,0x0BC4,\
0x0C8C,0x0D54,0x0E1C,0x0EE3,0x0FAB,0x1072,0x113A,0x1201,\
0x12C8,0x138F,0x1455,0x151C,0x15E2,0x16A8,0x176E,0x1833,\
0x18F9,0x19BE,0x1A82,0x1B47,0x1C0B,0x1CCF,0x1D93,0x1E57,\
0x1F1A,0x1FDD,0x209F,0x2161,0x2223,0x22E5,0x23A6,0x2467,\
0x2528,0x25E8,0x26A8,0x2767,0x2826,0x28E5,0x29A3,0x2A61,\
0x2B1F,0x2BDC,0x2C99,0x2D55,0x2E11,0x2ECC,0x2F87,0x3041,\
0x30FB,0x31B5,0x326E,0x3326,0x33DF,0x3496,0x354D,0x3604,\
0x36BA,0x376F,0x3824,0x38D9,0x398C,0x3A40,0x3AF2,0x3BA5,\
0x3C56,0x3D07,0x3DB8,0x3E68,0x3F17,0x3FC5,0x4073,0x4121,\
0x41CE,0x427A,0x4325,0x43D0,0x447A,0x4524,0x45CD,0x4675,\
0x471C,0x47C3,0x4869,0x490F,0x49B4,0x4A58,0x4AFB,0x4B9D,\
0x4C3F,0x4CE0,0x4D81,0x4E20,0x4EBF,0x4F5D,0x4FFB,0x5097,\
0x5133,0x51CE,0x5268,0x5302,0x539B,0x5432,0x54C9,0x5560,\
0x55F5,0x568A,0x571D,0x57B0,0x5842,0x58D3,0x5964,0x59F3,\
0x5A82,0x5B0F,0x5B9C,0x5C28,0x5CB3,0x5D3E,0x5DC7,0x5E4F,\
0x5ED7,0x5F5D,0x5FE3,0x6068,0x60EB,0x616E,0x61F0,0x6271,\
0x62F1,0x6370,0x63EE,0x646C,0x64E8,0x6563,0x65DD,0x6656,\
0x66CF,0x6746,0x67BC,0x6832,0x68A6,0x6919,0x698B,0x69FD,\
0x6A6D,0x6ADC,0x6B4A,0x6BB7,0x6C23,0x6C8E,0x6CF8,0x6D61,\
0x6DC9,0x6E30,0x6E96,0x6EFB,0x6F5E,0x6FC1,0x7022,0x7083,\
0x70E2,0x7140,0x719D,0x71F9,0x7254,0x72AE,0x7307,0x735E,\
0x73B5,0x740A,0x745F,0x74B2,0x7504,0x7555,0x75A5,0x75F3,\
0x7641,0x768D,0x76D8,0x7722,0x776B,0x77B3,0x77FA,0x783F,\
0x7884,0x78C7,0x7909,0x794A,0x7989,0x79C8,0x7A05,0x7A41,\
0x7A7C,0x7AB6,0x7AEE,0x7B26,0x7B5C,0x7B91,0x7BC5,0x7BF8,\
0x7C29,0x7C59,0x7C88,0x7CB6,0x7CE3,0x7D0E,0x7D39,0x7D62,\
0x7D89,0x7DB0,0x7DD5,0x7DFA,0x7E1D,0x7E3E,0x7E5F,0x7E7E,\
0x7E9C,0x7EB9,0x7ED5,0x7EEF,0x7F09,0x7F21,0x7F37,0x7F4D,\
0x7F61,0x7F74,0x7F86,0x7F97,0x7FA6,0x7FB4,0x7FC1,0x7FCD,\
0x7FD8,0x7FE1,0x7FE9,0x7FF0,0x7FF5,0x7FF9,0x7FFD,0x7FFE}

#define SIN_MASK        0x0300u
#define U0_90           0x0200u
#define U90_180         0x0300u
#define U180_270        0x0000u
#define U270_360        0x0100u

#define divSQRT_3 (int32_t)0x49E6    /* 1/sqrt(3) in q1.15 format=0.5773315*/

/*****************************************************************************
 * Module Preprocessor Macros
 ******************************************************************************/

/*****************************************************************************
 * Module Typedefs
 ******************************************************************************/

/*****************************************************************************
 * Function Prototypes
 ******************************************************************************/

/*****************************************************************************
 * Module Variable Definitions
 ******************************************************************************/

static const int16_t sin_cos_table[256] = SIN_COS_TABLE;

/*****************************************************************************
 * Function Definitions
 ******************************************************************************/

/**
 * @brief This function return the integer part of a fixed-point number
 *
 * @param f: fixed point number to convert
 * @param q: q format of the fixed point number as input
 * @return integer part of the fixed-point number
 */
int32_t lw_math_fix_2_int(fix_t f, uint32_t q) {

  if (f >= 0) {
      return (int32_t)(f / (1L << q));
  } else {
      return (int32_t)((f - (1L << q) + 1) / (1L << q));
  }
}

/**
 * @brief This function returns the integer part of a fixed-point
 *        number rounded of the nearest integer
 *
 * @param f: fixed point number to convert
 * @param q: q format of the fixed point number as input
 * @return integer part of the fixed-point number
 */
int32_t lw_math_fix_2_int_round(fix_t f, uint32_t q) {
  return lw_math_fix_2_int(f + (1L << q) / 2, q);
}

/**
 * @brief This function return the fractional part of a fixed-point number
 *
 * @param f: fixed point number from where to extract fractional part
 * @param q: q format of the fixed point number as input
 * @return fractional part of the fixed-point number
 */
fix_t lw_math_fix_fract_part(fix_t f, uint32_t q) {
  return f & ((1L << q) - 1);
}

/**
 * @brief  This function returns cosine and sine functions of the angle fed in
 *         input
 * @param  angle: angle in q1.15 format
 * @retval trig_components_t Cos(angle) and Sin(angle) in trig_components_t format
 */
trig_components_t lw_math_trig_functions(int16_t angle) {

  int32_t shindex;
  uint16_t uhindex;

  trig_components_t local_components;

  /* 10 bit index computation  */
  shindex = ((int32_t)32768 + (int32_t)angle);
  uhindex = (uint16_t)shindex;
  uhindex /= (uint16_t)64;

  switch((uint16_t)(uhindex) & SIN_MASK) {
    case U0_90: {
      local_components.sin = sin_cos_table[(uint8_t)(uhindex)];
      local_components.cos = sin_cos_table[(uint8_t)(0xFFu - (uint8_t)(uhindex))];
      break;
    }

    case U90_180: {
      local_components.sin = sin_cos_table[(uint8_t)(0xFFu - (uint8_t)(uhindex))];
      local_components.cos = -sin_cos_table[(uint8_t)(uhindex)];
      break;
    }

    case U180_270: {
      local_components.sin = -sin_cos_table[(uint8_t)(uhindex)];
      local_components.cos = -sin_cos_table[(uint8_t)(0xFFu - (uint8_t)(uhindex))];
      break;
    }

    case U270_360: {
      local_components.sin = -sin_cos_table[(uint8_t)(0xFFu - (uint8_t)(uhindex))];
      local_components.cos = sin_cos_table[(uint8_t)(uhindex)];
      break;
    }

    default: {
      break;
    }
  }

  return (local_components);
}

/**
 * @brief  It calculates the square root of a non-negative s32. It returns 0
 *         for negative s32.
 * @param  input int32_t number
 * @retval int32_t Square root of input (0 if input < 0)
 */
int32_t lw_math_sqrt(int32_t input) {

  uint8_t biter = 0u;
  int32_t wtemproot;
  int32_t wtemprootnew;

  if(input > 0) {

    if(input <= (int32_t)2097152) {
      wtemproot = (int32_t)128;
    }
    else {
      wtemproot = (int32_t)8192;
    }

    do {
      wtemprootnew = (wtemproot + input / wtemproot) / (int32_t)2;
      if(wtemprootnew == wtemproot) {
        biter = 6u;
      }
      else {
        biter++;
        wtemproot = wtemprootnew;
      }
    }while (biter < 6u);
  }
  else {
    wtemprootnew = (int32_t)0;
  }

  return (wtemprootnew);
}

/**
  * @brief  This function transforms components a and b (which are
  *         directed along axes each displaced by 120 degrees) into components
  *         alpha and beta.
  *                               alpha = Ia
  *                       beta = -(2 * Ib + Ia) / sqrt(3)
  * @param  input: component a and b in ab_t format
  * @retval Components alpha and beta in alphabeta_t format
  */
alphabeta_t lw_math_clarke(ab_t input) {

  alphabeta_t output;

  int32_t a_divSQRT3_tmp;
  int32_t b_divSQRT3_tmp;
  int32_t wbeta_tmp;
  int16_t hbeta_tmp;

  /* qIalpha = qIas*/
  output.alpha = input.a;

  a_divSQRT3_tmp = divSQRT_3 * ((int32_t)input.a);

  b_divSQRT3_tmp = divSQRT_3 * ((int32_t)input.b);

  /*qIbeta = -(2*qIbs+qIas)/sqrt(3)*/

  /* WARNING: the below instruction is not MISRA compliant, user should verify
    that Cortex-M3 assembly instruction ASR (arithmetic shift right) is used by
    the compiler to perform the shift (instead of LSR logical shift right) */
  //cstat !MISRAC2012-Rule-1.3_n !ATH-shift-neg !MISRAC2012-Rule-10.1_R6
  /*wbeta_tmp = (-(a_divSQRT3_tmp) - (b_divSQRT3_tmp) - (b_divSQRT3_tmp)) >> 15;*/

  wbeta_tmp = (-(a_divSQRT3_tmp) - (b_divSQRT3_tmp) - (b_divSQRT3_tmp)) / 32768;


  /* Check saturation of Ibeta */
  if (wbeta_tmp > INT16_MAX)
  {
    hbeta_tmp = INT16_MAX;
  }
  else if (wbeta_tmp < (-32768))
  {
    hbeta_tmp =  ((int16_t)-32768);
  }
  else
  {
    hbeta_tmp = ((int16_t)wbeta_tmp);
  }

  output.beta = hbeta_tmp;

  if (((int16_t )-32768) == output.beta)
  {
    output.beta = -32767;
  }

  return (output);
}

/**
  * @brief  This function transforms components alpha and beta, which
  *         belong to a stationary qd reference frame, to a rotor flux
  *         synchronous reference frame (properly oriented), so as q and d.
  *                   d= alpha *sin(theta) + beta * cos(theta)
  *                   q= alpha *cos(theta) - beta * sin(theta)
  * @param  input: components values alpha and beta in alphabeta_t format
  * @param  theta: rotating frame angular position in q1.15 format
  * @retval Components q and d in qd_t format
  */
qd_t lw_math_park(alphabeta_t input, int16_t theta) {

  qd_t output;
  int32_t d_tmp_1;
  int32_t d_tmp_2;
  int32_t q_tmp_1;
  int32_t q_tmp_2;
  int32_t wqd_tmp;
  int16_t hqd_tmp;
  trig_components_t Local_Vector_Components;

  Local_Vector_Components = lw_math_trig_functions(theta);

  /*No overflow guaranteed*/
  q_tmp_1 = input.alpha * ((int32_t )Local_Vector_Components.cos);

  /*No overflow guaranteed*/
  q_tmp_2 = input.beta * ((int32_t)Local_Vector_Components.sin);

  /*Iq component in Q1.15 Format */
  /* WARNING: the below instruction is not MISRA compliant, user should verify
    that Cortex-M3 assembly instruction ASR (arithmetic shift right) is used by
    the compiler to perform the shift (instead of LSR logical shift right) */
  //cstat !MISRAC2012-Rule-1.3_n !ATH-shift-neg !MISRAC2012-Rule-10.1_R6
  /*wqd_tmp = (q_tmp_1 - q_tmp_2) >> 15;  */

  wqd_tmp = (q_tmp_1 - q_tmp_2) / 32768;

  /* Check saturation of Iq */
  if (wqd_tmp > INT16_MAX)
  {
    hqd_tmp = INT16_MAX;
  }
  else if (wqd_tmp < (-32768))
  {
    hqd_tmp = ((int16_t)-32768);
  }
  else
  {
    hqd_tmp = ((int16_t)wqd_tmp);
  }

  output.q = hqd_tmp;

  if (((int16_t )-32768) == output.q)
  {
    output.q = -32767;
  }

  /*No overflow guaranteed*/
  d_tmp_1 = input.alpha * ((int32_t )Local_Vector_Components.sin);

  /*No overflow guaranteed*/
  d_tmp_2 = input.beta * ((int32_t )Local_Vector_Components.cos);

  /*Id component in Q1.15 Format */
  /* WARNING: the below instruction is not MISRA compliant, user should verify
    that Cortex-M3 assembly instruction ASR (arithmetic shift right) is used by
    the compiler to perform the shift (instead of LSR logical shift right) */
  //cstat !MISRAC2012-Rule-1.3_n !ATH-shift-neg !MISRAC2012-Rule-10.1_R6
  /* wqd_tmp = (d_tmp_1 + d_tmp_2) >> 15; */

  wqd_tmp = (d_tmp_1 + d_tmp_2) / 32768;

  /* Check saturation of Id */
  if (wqd_tmp > INT16_MAX)
  {
    hqd_tmp = INT16_MAX;
  }
  else if (wqd_tmp < (-32768))
  {
    hqd_tmp = ((int16_t)-32768);
  }
  else
  {
    hqd_tmp = ((int16_t)wqd_tmp);
  }

  output.d = hqd_tmp;

  if (((int16_t)-32768) == output.d)
  {
    output.d = -32767;
  }

  return (output);
}

/**
  * @brief  This function transforms the input component q and d, to a stationary reference
  *         frame, so as to obtain alpha and beta:
  *                  alfa= q * cos(theta)+ d * sin(theta)
  *                  beta= -q * sin(theta)+ d * cos(theta)
  * @param  input: input component q and d in qd_t format
  * @param  theta: angular position in q1.15 format
  * @retval output component alpha and beta in alphabeta_t format
  */
alphabeta_t lw_math_rev_park(qd_t input, int16_t theta) {

  int32_t alpha_tmp1;
  int32_t alpha_tmp2;
  int32_t beta_tmp1;
  int32_t beta_tmp2;
  trig_components_t Local_Vector_Components;
  alphabeta_t output;

  Local_Vector_Components = lw_math_trig_functions(theta);

  /*No overflow guaranteed*/
  alpha_tmp1 = input.q * ((int32_t)Local_Vector_Components.cos);
  alpha_tmp2 = input.d * ((int32_t)Local_Vector_Components.sin);


  /* WARNING: the below instruction is not MISRA compliant, user should verify
    that Cortex-M3 assembly instruction ASR (arithmetic shift right) is used by
    the compiler to perform the shift (instead of LSR logical shift right) */

  //cstat !MISRAC2012-Rule-1.3_n !ATH-shift-neg !MISRAC2012-Rule-10.1_R6
  /*output.alpha = (int16_t)(((alpha_tmp1) + (alpha_tmp2)) >> 15);*/

  output.alpha = (int16_t)(((alpha_tmp1) + (alpha_tmp2)) / 32768);


  beta_tmp1 = input.q * ((int32_t)Local_Vector_Components.sin);
  beta_tmp2 = input.d * ((int32_t)Local_Vector_Components.cos);

  /* WARNING: the below instruction is not MISRA compliant, user should verify
  that Cortex-M3 assembly instruction ASR (arithmetic shift right) is used by
  the compiler to perform the shift (instead of LSR logical shift right) */
  //cstat !MISRAC2012-Rule-1.3_n !ATH-shift-neg !MISRAC2012-Rule-10.1_R6
  /*output.beta = (int16_t)((beta_tmp2 - beta_tmp1) >> 15);*/

  output.beta = (int16_t)((beta_tmp2 - beta_tmp1) / 32768);


  return (output);
}

/*************** END OF FUNCTIONS ********************************************/


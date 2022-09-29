/*****************************************************************************
 * Filename              :   lw_math.h
 * Author                :   Giulio Dalla Vecchia
 * Origin Date           :   5 may 2022
 *
 * Copyright (c) 2022 Giulio Dalla Vecchia. All rights reserved.
 *
 ******************************************************************************/

/** @file lw_math.h
 *  @brief This module declares an interface to control the basic math
 *         operations in fixed-point format
 */

#ifndef LW_MATH_H_
#define LW_MATH_H_

/*****************************************************************************
 * Includes
 ******************************************************************************/
#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

/**
 * \defgroup        lw_math
 * \brief           Lightweight Mathematical Library
 * \{
 */

/*****************************************************************************
 * Module Preprocessor Constants
 ******************************************************************************/

/* convert to and from integer */
#define INT2FIX(d, q) ((fix_t)((d) << (q)))
#define FIX2INT(d, q) ((int32_t)((d) >> (q)))

/* convert to and from floating point */
#define FLT2FIX(d, q) ((fix_t)((double)(d) * (1<<(q))))
#define FIX2FLT(a, q) ((double)(a) / (double)(1<<(q)))

/* The basic operations perfomed on two numbers a and b of fixed
 point q format returning the answer in q format */
#define FADD(a,b) ((a)+(b))
#define FSUB(a,b) ((a)-(b))
#define FMUL(a,b,q) ((fix_t)(((int64_t)(a)*(int64_t)(b))>>(q)))
#define FDIV(a,b,q) ((fix_t)(((int64_t)(a)<<(q))/(int64_t)(b)))

/* The basic operations where a is of fixed point q format and b is
 an integer */
#define FADDI(a,b,q) ((a)+((b)<<(q)))
#define FSUBI(a,b,q) ((a)-((b)<<(q)))
#define FMULI(a,b) ((a)*(b))
#define FDIVI(a,b) ((a)/(b))

/* convert a from q1 format to q2 format */
#define FCONV(a, q1, q2) (((q2)>(q1)) ? (a)<<((q2)-(q1)) : (a)>>((q1)-(q2)))

/* the general operation between a in q1 format and b in q2 format
 returning the result in q3 format */
#define FADDG(a,b,q1,q2,q3) (FCONV(a,q1,q3)+FCONV(b,q2,q3))
#define FSUBG(a,b,q1,q2,q3) (FCONV(a,q1,q3)-FCONV(b,q2,q3))
#define FMULG(a,b,q1,q2,q3) FCONV((a)*(b), (q1)+(q2), q3)
#define FDIVG(a,b,q1,q2,q3) (FCONV(a, q1, (q2)+(q3))/(b))

/* Square root of a number in q format */
#define FSQRT(a, q) (lw_math_sqrt(a << q))

/*****************************************************************************
 * Module Preprocessor Macros
 ******************************************************************************/

/*****************************************************************************
 * Module Typedefs
 ******************************************************************************/

/**< Fixed point type (as example could be s16.15) */
typedef int32_t fix_t;

/**
 * @brief  Trigonometrical functions type definition
 */
typedef struct {
  int16_t cos;
  int16_t sin;
} trig_components_t;

/**
 * @brief Two components a,b type definition
 */
typedef struct {
  int16_t a;
  int16_t b;
} ab_t;

/**
 * @brief Two components q, d type definition
 */
typedef struct {
  int16_t q;
  int16_t d;
} qd_t;

/**
 * @brief Two components alpha, beta type definition
 */
typedef struct {
  int16_t alpha;
  int16_t beta;
} alphabeta_t;

/*****************************************************************************
 * Module Variable Definitions
 ******************************************************************************/

/*****************************************************************************
 * Function Prototypes
 ******************************************************************************/

/**
 * @brief This function return the integer part of a fixed-point number
 *
 * @param f: fixed point number to convert
 * @param q: q format of the fixed point number as input
 * @return integer part of the fixed-point number
 */
int32_t lw_math_fix_2_int(fix_t f, uint32_t q);

/**
 * @brief This function returns the integer part of a fixed-point
 *        number rounded of the nearest integer
 *
 * @param f: fixed point number to convert
 * @param q: q format of the fixed point number as input
 * @return integer part of the fixed-point number
 */
int32_t lw_math_fix_2_int_round(fix_t f, uint32_t q);

/**
 * @brief This function return the fractional part of a fixed-point number
 *
 * @param f: fixed point number from where to extract fractional part
 * @param q: q format of the fixed point number as input
 * @return fractional part of the fixed-point number
 */
fix_t lw_math_fix_fract_part(fix_t f, uint32_t q);

/**
 * @brief  This function returns cosine and sine functions of the angle fed in
 *         input
 * @param  angle: angle in q1.15 format
 * @retval trig_components_t Cos(angle) and Sin(angle) in trig_components_t format
 */
trig_components_t lw_math_trig_functions(int16_t angle);

/**
 * @brief  It calculates the square root of a non-negative s32. It returns 0
 *         for negative s32.
 * @param  input int32_t number
 * @retval int32_t Square root of input (0 if input < 0)
 */
int32_t lw_math_sqrt(int32_t input);

/**
 * @brief  This function transforms components a and b (which are
 *         directed along axes each displaced by 120 degrees) into components
 *         alpha and beta.
 *                               alpha = Ia
 *                       beta = -(2 * Ib + Ia) / sqrt(3)
 * @param  input: component a and b in ab_t format
 * @retval Components alpha and beta in alphabeta_t format
 */
alphabeta_t lw_math_clarke(ab_t input);

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
qd_t lw_math_park(alphabeta_t input, int16_t theta);

/**
 * @brief  This function transforms the input component q and d, to a stationary reference
 *         frame, so as to obtain alpha and beta:
 *                  alfa= q * cos(theta)+ d * sin(theta)
 *                  beta= -q * sin(theta)+ d * cos(theta)
 * @param  input: input component q and d in qd_t format
 * @param  theta: angular position in q1.15 format
 * @retval output component alpha and beta in alphabeta_t format
 */
alphabeta_t lw_math_rev_park(qd_t input, int16_t theta);

/**
 * \}
 */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /*LW_MATH_H_*/

/*** End of File *************************************************************/

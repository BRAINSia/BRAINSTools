#if 0
/*******************************************************************************
  Common point for the registration parameters

  Sequence of parameters for affine:
  Translation x, y, z
  Rotation x, y, z
  Scaling x, y, z
  Skew x, y, z
 ***************************************************************************** */

#ifndef _RegistrationParameters_h
#define _RegistrationParameters_h

#define MU_REGISTER_MAJORVER "2"
#define MU_REGISTER_MINORVER "0a"

//
// Line search steps
//

#define MU_AFFINE_STEP_TRANSLATE 10.0
#define MU_AFFINE_STEP_ROTATE 0.1
#define MU_AFFINE_STEP_SCALE 0.1
#define MU_AFFINE_STEP_SKEW 0.01

//
// Optimization order
//

// Translation, rotation, scale, then skew

#define MU_AFFINE_ORDER0 0
#define MU_AFFINE_ORDER1 1
#define MU_AFFINE_ORDER2 2
#define MU_AFFINE_ORDER3 3
#define MU_AFFINE_ORDER4 4
#define MU_AFFINE_ORDER5 5
#define MU_AFFINE_ORDER6 6
#define MU_AFFINE_ORDER7 7
#define MU_AFFINE_ORDER8 8
#define MU_AFFINE_ORDER9 9
#define MU_AFFINE_ORDER10 10
#define MU_AFFINE_ORDER11 11

/*
// Translation, scale, rotation, then skew

#define MU_AFFINE_ORDER0 0
#define MU_AFFINE_ORDER1 1
#define MU_AFFINE_ORDER2 2
#define MU_AFFINE_ORDER3 6
#define MU_AFFINE_ORDER4 7
#define MU_AFFINE_ORDER5 8
#define MU_AFFINE_ORDER6 3
#define MU_AFFINE_ORDER7 4
#define MU_AFFINE_ORDER8 5
#define MU_AFFINE_ORDER9 9
#define MU_AFFINE_ORDER10 10
#define MU_AFFINE_ORDER11 11
*/

#endif
#endif

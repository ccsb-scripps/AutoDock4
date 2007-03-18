#ifndef _REP_CONSTANTS
#define _REP_CONSTANTS

// Translation
#define X_TRANSLATION_INDEX 0
#define Z_TRANSLATION_INDEX 2

// Quaternion
#define QX_ROTATION_INDEX 3
#define QY_ROTATION_INDEX 4
#define QZ_ROTATION_INDEX 5
#define QW_ROTATION_INDEX 6

// Axis-Angle
#define X_ROTATION_INDEX 3
#define Y_ROTATION_INDEX 4
#define Z_ROTATION_INDEX 5

#define ROTATION_ANGLE_INDEX QW_ROTATION_INDEX

#define is_translation_index(i) (((i) >= X_TRANSLATION_INDEX) && ((i) <= Z_TRANSLATION_INDEX))
#define is_axis_index(i) (((i) >= X_ROTATION_INDEX) && ((i) <= Z_ROTATION_INDEX))
#define is_angle_index(i) ((i) == ROTATION_ANGLE_INDEX)
#define is_rotation_index(i) (((i) >= QX_ROTATION_INDEX) && ((i) <= QW_ROTATION_INDEX))

#endif

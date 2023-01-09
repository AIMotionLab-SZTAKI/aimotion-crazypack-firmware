/*
The MIT License (MIT)

Copyright (c) 2018 Wolfgang Hoenig and James Alan Preiss

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <math.h>

#include "param.h"
#include "log.h"
#include "math3d.h"
#include "position_controller.h"
#include "controller_pid.h"
#include "controller_flip.h"
#include "controller_geom.h"
#include "stabilizer.h"

#define GRAVITY_MAGNITUDE (9.81f)

static float g_vehicleMass = 0.032; // TODO: should be CF global for other modules
// static float massThrust = 132000;

// Inertia matrix components
static float Ixx = 0.0000158;
// static float Iyy = 0.0000158;
// static float Izz = 0.00002926;

// Flip parameters
static float U1 = 0.5259;
static float Th1 = -42.346;
static float T1 = 0.08219;

static float U2 = 0.3795;
static float Th2 = 297.346;
static float T2 = 0.22265;

static float U3 = 0.17489;
static float Th3 = 0.0;
static float T3 = 0.1209;

static float U4 = 0.3795;
static float Th4 = -297.346;
static float T4 = 0.22193;

static float U5 = 0.50265;
static float Th5 = 59.2655;
static float T5 = 0.05508;

static float t = 0;

static float cmd_thrust;
static float cmd_roll;
static float cmd_pitch;
static float cmd_yaw;
static float r_roll;
static float r_pitch;
static float r_yaw;
static float accelz;

static bool isFlipControl = false;
static bool wasFlipControl = false;

static int controller_type = 1; // 1: PID, 2: Mellinger

static float x_bound = 1;
static float y_bound = 1;
static float z_bound = 0.2;
static float vel_bound = 0.5;


void controllerFlipFFReset(void)
{
    t = 0;
}

void controllerFlipFFInit(void)
{
  controllerFlipFFReset();
  controllerPidInit();          // we use the PID controller
  controllerGeomInit();
  t = 0;
}

bool controllerFlipFFTest(void)
{
  return true;
}

void controllerFlipFF(control_t *control, setpoint_t *setpoint,
                                         const sensorData_t *sensors,
                                         const state_t *state,
                                         const uint32_t tick)
{
  float dt;

  if (!RATE_DO_EXECUTE(ATTITUDE_RATE, tick)) {
    return;
  }

  dt = (float)(1.0f/ATTITUDE_RATE);

  if (!isFlipControl) {
    if (wasFlipControl)
    {
      if(((float)fabs(state->position.x) > x_bound) || ((float)fabs(state->position.y) > y_bound || state->position.z < z_bound))
      {
        stabilizerSetEmergencyStop();
      }
      if((float)fabs(state->velocity.z) < vel_bound)
      {
        //wasFlipControl = false;
      }
      /*setpoint->position.z = recovery_setpoint;
      setpoint->mode.x = modeDisable;
      setpoint->mode.y = modeDisable;
      setpoint->mode.yaw = modeDisable;
      setpoint->attitude.roll = 0;
      setpoint->attitude.pitch = 0;*/
    }
    switch (controller_type)
    {
    case 1:
      controllerGeom(control, setpoint, sensors, state, tick);
      break;
    case 2:
      controllerPid(control, setpoint, sensors, state, tick);
    default:
      break;
    }
  } else {
    if (t > T1 + T2 + T3 + T4 + T5) {
    isFlipControl = false;
    wasFlipControl = true;
    } else if (t > T1 + T2 + T3 + T4) {
    control->thrust = U5 * g_vehicleMass;
    control->pitch = Ixx * Th5;
    } else if (t > T1 + T2 + T3) {
    control->thrust = U4 * g_vehicleMass;
    control->pitch = Ixx * Th4;
    } else if (t > T1 + T2) {
    control->thrust = U3 * g_vehicleMass;
    control->pitch = Ixx * Th3;
    } else if (t > T1) {
    control->thrust = U2 * g_vehicleMass;
    control->pitch = Ixx * Th2;
    } else {
    control->thrust = U1 * g_vehicleMass;
    control->pitch = Ixx * Th1;
    }
    t += dt;

    cmd_thrust = control->thrust;
    r_roll = radians(sensors->gyro.x);
    r_pitch = -radians(sensors->gyro.y);
    r_yaw = radians(sensors->gyro.z);
    accelz = sensors->acc.z;

    if (control->thrust > 0) {                                 
        control->roll = 0;
        //control->pitch = M.x*(float)(4.7784e+5 / (2.0 * 0.046)) + 0*M.y;
        control->yaw = 0;
        // control->roll = 1000; // clamp(M.x, -32000, 32000);
        // control->pitch = 0; //clamp(M.y, -32000, 32000);
        // control->yaw = 0; // clamp(-M.z, -32000, 32000);

    } else {
        control->roll = 0;
        control->pitch = 0;
        control->yaw = 0;
        controllerFlipFFReset();
    }
  }
  cmd_roll = control->roll;
  cmd_pitch = control->pitch;
  cmd_yaw = control->yaw;
  cmd_thrust = control->thrust;
}

PARAM_GROUP_START(ctrlFlip)
PARAM_ADD(PARAM_FLOAT, U1, &U1)
PARAM_ADD(PARAM_FLOAT, Th1, &Th1)
PARAM_ADD(PARAM_FLOAT, T1, &T1)
PARAM_ADD(PARAM_FLOAT, U2, &U2)
PARAM_ADD(PARAM_FLOAT, Th2, &Th2)
PARAM_ADD(PARAM_FLOAT, T2, &T2)
PARAM_ADD(PARAM_FLOAT, U3, &U3)
PARAM_ADD(PARAM_FLOAT, Th3, &Th3)
PARAM_ADD(PARAM_FLOAT, T3, &T3)
PARAM_ADD(PARAM_FLOAT, U4, &U4)
PARAM_ADD(PARAM_FLOAT, Th4, &Th4)
PARAM_ADD(PARAM_FLOAT, T4, &T4)
PARAM_ADD(PARAM_FLOAT, U5, &U5)
PARAM_ADD(PARAM_FLOAT, Th5, &Th5)
PARAM_ADD(PARAM_FLOAT, T5, &T5)
PARAM_GROUP_STOP(ctrlFlip)

LOG_GROUP_START(ctrlFlip)
LOG_ADD(LOG_FLOAT, cmd_thrust, &cmd_thrust)
LOG_ADD(LOG_FLOAT, cmd_roll, &cmd_roll)
LOG_ADD(LOG_FLOAT, cmd_pitch, &cmd_pitch)
LOG_ADD(LOG_FLOAT, cmd_yaw, &cmd_yaw)
LOG_ADD(LOG_FLOAT, r_roll, &r_roll)
LOG_ADD(LOG_FLOAT, r_pitch, &r_pitch)
LOG_ADD(LOG_FLOAT, r_yaw, &r_yaw)
LOG_ADD(LOG_FLOAT, accelz, &accelz)
LOG_ADD(LOG_FLOAT, t, &t)
LOG_GROUP_STOP(ctrlFlip)
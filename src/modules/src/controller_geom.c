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

/*
Geometric tracking controller.
*/

#include <math.h>

#include "param.h"
#include "log.h"
#include "math3d.h"
#include "stabilizer.h"
#include "physicalConstants.h"
#include "controller_geom.h"
#include "controller_flip.h"
#include "interpolation.h"


// Inertia matrix components
static float Ixx = 0.000014;
static float Izz = 0.0000217;

static float t = 0;
static float dt = (float)(1.0f/ATTITUDE_RATE);

const static float thrust_scale = 132000; // 119460.0; TODO: find a solution (scale gains somehow)
const static float l = 0.0325; // L/sqrt(2), where L is half prop-to-prop length
const static float b = 0.025; // b/k to scale yaw

static float kr_xy = 0.4;
static float kr_z = 1.25;
static float kv_xy = 0.2;
static float kv_z = 0.4;
static float kR_xy = 0.017; // 70000/thrust_scale*l instead of 70000/rollpitch_scale;
static float kR_z = 0.011; // 60000/thrust_scale*b instead of 60000/yaw_scale
static float kw_xy = 0.0049; // 20000/thrust_scale*l instead of 20000/rollpitch_scale;
static float kw_z = 0.00225; // 12000/thrust_scale*b instead of 12000/yaw_scale
static float kd_omega_rp = 0.000049; // 200/thrust_scale*l instead of 200/rollpitch_scale

// Logging variables

static float cmd_thrust_g;
static float cmd_roll;
static float cmd_pitch;
static float cmd_yaw;
static float r_roll;
static float r_pitch;
static float r_yaw;

static float q0;
static float q2;
static struct quat q;
static float timescale = 1;
static struct vec eR, ew, wd;
static float p_des, p;

static float prev_omega_roll;
static float prev_omega_pitch;
static float prev_setpoint_omega_roll;
static float prev_setpoint_omega_pitch;
static float err_d_roll = 0;
static float err_d_pitch = 0;

static bool mode = false;  // 0: position control, 1: attitude control
static float thrust_tmp;

static float g_vehicleMass = 0.027;

static float psi = 0;

//static bool gp = false;
//static float u0 = 0;
//static float u1 = 0;
//static float gr0[5] = {-0.05 , -0.025,  0.   ,  0.025,  0.05};
//static float gr1[9] = {-1.  , -0.85, -0.7 , -0.55, -0.4 , -0.25, -0.1 ,  0.05,  0.2};
//static float gr2[5] = {-0.2, -0.1,  0. ,  0.1,  0.2};
//static float gr3[9] = {-0.2   ,  0.1125,  0.425 ,  0.7375,  1.05  ,  1.3625,  1.675 ,  1.9875,  2.3};
//static int min0_i = 0, min1_i = 0, min2_i = 0, min3_i = 0;

/*static float u2 = 0;
static float sf1 = 0.135f;
static float sf2 = 0.477f;
static float X[48][4] = {{0.000,-0.002,-0.002,-0.003},
{0.000,-0.002,0.013,0.040},
{0.003,-0.012,0.043,0.142},
{0.010,-0.035,0.036,0.171},
{0.014,-0.060,0.009,0.164},
{0.016,-0.085,0.007,0.183},
{0.017,-0.115,0.011,0.251},
{0.018,-0.156,0.008,0.341},
{0.019,-0.211,0.007,0.465},
{0.018,-0.317,-0.058,0.737},
{0.011,-0.428,-0.042,0.937},
{0.009,-0.559,-0.015,1.191},
{0.009,-0.707,-0.023,1.496},
{0.009,-0.854,-0.057,1.852},
{0.012,-0.970,-0.056,2.125},
{0.016,-0.997,-0.013,2.144},
{0.020,-0.925,0.018,1.968},
{0.022,-0.782,0.010,1.680},
{0.025,-0.614,-0.003,1.402},
{0.028,-0.446,-0.009,1.156},
{0.028,-0.298,0.062,0.881},
{0.016,-0.182,0.105,0.581},
{0.009,-0.108,-0.074,0.298},
{0.018,-0.077,0.002,0.038},
{0.000,-0.000,0.000,0.001},
{0.001,-0.002,0.017,0.039},
{0.006,-0.014,0.053,0.134},
{0.016,-0.043,0.029,0.149},
{0.019,-0.065,0.023,0.150},
{0.022,-0.091,0.027,0.203},
{0.026,-0.127,0.025,0.272},
{0.028,-0.172,0.010,0.352},
{0.028,-0.234,0.003,0.490},
{0.024,-0.321,-0.074,0.747},
{0.014,-0.437,-0.050,0.916},
{0.012,-0.568,0.001,1.135},
{0.013,-0.666,-0.007,1.473},
{0.011,-0.865,-0.024,1.818},
{0.010,-0.974,-0.013,2.103},
{0.010,-0.996,0.018,2.150},
{0.007,-0.920,0.018,2.000},
{0.006,-0.829,0.006,1.706},
{0.004,-0.610,-0.005,1.406},
{0.002,-0.446,-0.005,1.185},
{-0.003,-0.295,0.038,0.990},
{-0.007,-0.157,-0.007,0.832},
{-0.009,-0.084,0.058,0.431},
{-0.023,-0.030,0.073,-0.027}};
static float alpha1[48] = {-1248.766,459.255,44.214,-34.284,127.946,-216.751,-2.501,121.021,-55.371,1.990,1.424,-9.429,-8.707,2.210,-1.681,2.783,-3.673,0.432,-1.284,2.732,-2.483,-1.006,4.522,27.951,1387.830,-645.825,17.690,14.348,39.608,-64.078,32.931,-17.654,7.302,3.362,-1.635,4.712,5.039,4.531,-2.717,-2.625,3.136,-2.654,6.305,3.962,-4.778,2.602,-6.563,-2.166};
static float alpha2[48] = {-106.574,176.405,-227.534,330.329,-243.287,289.092,-259.870,220.247,-88.794,17.381,1.528,2.153,5.764,-5.659,6.487,-1.496,-0.003,-0.692,-2.744,0.853,-2.063,-0.351,-4.260,11.433,-73.246,108.075,3.004,-354.976,362.757,-250.661,161.215,-106.161,45.830,-10.685,-5.976,4.136,-2.143,6.399,-21.441,10.587,7.035,-3.006,-1.563,-6.298,3.222,-4.217,7.658,0.271};
static float lam1[4] = {16119.646,11.104,222.851,9.934};
static float lam2[4] = {7222.721,1.805,54.594,13.825};*/

//static int ltb0[2025] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21,-54,-174,-115,17,0,0,0,0,21,-54,-175,-115,17,0,0,0,0,21,-54,-175,-115,17,0,0,0,0,21,-54,-174,-114,17,0,0,0,0,21,-54,-173,-113,16,1,-5,0,-2,86,-109,-177,64,22,1,-5,0,-2,86,-110,-178,64,22,1,-5,0,-2,87,-110,-178,64,22,1,-5,0,-2,86,-109,-178,64,22,1,-5,0,-2,86,-109,-176,63,22,21,-79,11,-11,216,-91,-129,150,15,21,-80,11,-11,217,-91,-130,151,15,21,-81,11,-11,218,-91,-130,151,15,21,-81,11,-11,217,-91,-130,150,15,21,-80,11,-11,215,-90,-129,149,15,197,-682,95,-25,334,10,-87,93,6,199,-689,96,-26,336,10,-88,94,6,199,-692,97,-26,336,10,-88,94,6,199,-693,97,-26,335,10,-88,93,6,198,-690,96,-25,332,10,-88,92,6,1070,-3422,485,-34,309,66,-48,28,1,1077,-3456,489,-34,311,66,-48,28,1,1080,-3474,491,-34,311,66,-48,28,1,1078,-3477,491,-34,310,66,-48,27,1,1071,-3465,488,-34,307,66,-48,27,1,3405,-10073,1446,-26,167,43,-16,4,0,3429,-10173,1458,-26,168,43,-16,4,0,3438,-10228,1464,-26,168,43,-16,4,0,3432,-10238,1464,-26,167,43,-16,4,0,3411,-10204,1457,-26,166,43,-16,4,0,6381,-17391,2537,-11,52,12,-3,0,0,6426,-17564,2559,-11,52,12,-3,0,0,6443,-17663,2569,-11,52,12,-3,0,0,6432,-17684,2569,-11,52,12,-3,0,0,6393,-17628,2557,-11,51,12,-3,0,0,7041,-17602,2619,-3,9,2,0,0,0,7090,-17780,2641,-3,9,2,0,0,0,7109,-17884,2652,-3,9,2,0,0,0,7097,-17907,2652,-3,9,2,0,0,0,7054,-17855,2640,-3,9,2,0,0,0,4574,-10437,1591,0,0,0,0,0,0,4606,-10545,1604,0,0,0,0,0,0,4619,-10608,1611,0,0,0,0,0,0,4611,-10624,1611,0,0,0,0,0,0,4583,-10594,1603,0,0,0,0,0,0,0,0,-3,-24,7,-859,-1126,-3684,432,0,0,-3,-24,7,-865,-1135,-3714,435,0,0,-3,-24,7,-867,-1140,-3728,436,0,0,-3,-24,7,-865,-1139,-3725,436,0,0,-3,-24,7,-860,-1134,-3707,433,1,-6,-66,-223,85,-3010,-2507,-6197,670,1,-6,-67,-223,86,-3029,-2529,-6248,675,1,-6,-67,-223,87,-3036,-2540,-6271,676,1,-6,-67,-222,87,-3029,-2540,-6268,675,1,-6,-66,-219,87,-3009,-2528,-6237,670,22,-95,-584,-1116,504,-6608,-2041,-6692,481,23,-95,-587,-1119,509,-6651,-2061,-6746,484,23,-96,-588,-1117,511,-6664,-2071,-6771,485,23,-96,-586,-1110,512,-6648,-2073,-6767,484,23,-95,-581,-1098,510,-6604,-2066,-6733,481,198,-820,-2581,-2814,1410,-9505,-271,-4274,183,202,-826,-2596,-2816,1424,-9565,-276,-4309,184,204,-829,-2600,-2806,1431,-9583,-280,-4325,184,206,-828,-2593,-2783,1431,-9559,-283,-4322,183,207,-824,-2574,-2748,1426,-9494,-285,-4300,182,1022,-4129,-5761,-2245,1445,-9304,433,-1574,38,1040,-4161,-5798,-2223,1459,-9363,435,-1587,38,1054,-4175,-5810,-2191,1467,-9380,435,-1592,38,1063,-4171,-5796,-2149,1468,-9356,434,-1591,38,1068,-4148,-5758,-2099,1463,-9293,430,-1583,38,3099,-12192,-5744,4339,-950,-6240,224,-332,4,3154,-12286,-5789,4418,-956,-6279,226,-335,4,3196,-12327,-5808,4479,-958,-6291,227,-336,4,3224,-12315,-5802,4520,-955,-6276,226,-336,4,3238,-12249,-5770,4540,-949,-6234,225,-334,4,5533,-21154,-622,11545,-3310,-2784,45,-40,0,5630,-21318,-642,11666,-3336,-2802,46,-40,0,5704,-21390,-659,11736,-3347,-2808,46,-40,0,5754,-21368,-674,11754,-3344,-2801,46,-40,0,5779,-21253,-686,11722,-3326,-2783,45,-40,0,5815,-21591,3328,10863,-2816,-793,5,-2,0,5917,-21758,3333,10960,-2839,-798,5,-2,0,5995,-21832,3324,11009,-2849,-800,5,-2,0,6047,-21810,3300,11010,-2847,-798,5,-2,0,6073,-21693,3262,10963,-2832,-793,5,-2,0,3598,-12969,2714,5211,-1139,-139,0,0,0,3662,-13070,2723,5254,-1148,-140,0,0,0,3710,-13115,2720,5275,-1152,-140,0,0,0,3742,-13102,2705,5273,-1151,-140,0,0,0,3758,-13032,2678,5247,-1145,-139,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//static int ltb1[2025] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-28,627,-1848,-745,0,0,0,0,0,-90,-542,-2781,-701,0,0,0,0,0,0,-10,-196,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-459,1134,-2482,-569,0,0,0,0,0,-1929,-6204,-3074,-1397,0,0,0,0,0,0,-58,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,-15,-1209,-1519,-386,-37,0,0,0,-13,-352,-8732,-12533,-867,-84,0,0,0,0,0,-7,-50,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,-335,-298,-619,-434,-17,0,0,0,68,-104,-6363,-12577,-5765,-117,0,0,0,0,0,-4,-57,-20,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-110,-1877,-356,-56,-10,0,0,-4,-44,-135,-1531,-8985,-8672,-950,-1,0,-2,3,-7,20,-7,-78,-4,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-778,-1602,-65,0,0,0,0,-520,-6436,-5313,-3546,-7958,-1540,-36,0,0,-274,88,-25,-181,-263,-30,0,0,0,-202,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,157,-2098,-235,-3,0,0,0,0,-6516,-21506,-3248,-3526,-3282,-67,0,0,0,-2891,-701,168,-772,-159,-2,0,0,0,-2135,-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,239,-456,-17,0,0,0,0,0,-9534,29071,62,-933,-175,-1,0,0,0,-3512,-2604,147,-220,-7,0,0,0,0,-2588,-62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,46,-1,0,0,0,0,0,-1622,10365,75,-35,-1,0,0,0,0,-491,-574,7,-6,0,0,0,0,0,-361,-8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-28,628,-1849,-745,0,0,0,0,0,-90,-542,-2782,-701,0,0,0,0,0,0,-10,-196,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-459,1134,-2482,-569,0,0,0,0,0,-1929,-6205,-3075,-1397,0,0,0,0,0,0,-58,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,-15,-1209,-1520,-386,-37,0,0,0,-13,-352,-8734,-12536,-867,-84,0,0,0,0,0,-7,-50,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,-335,-298,-619,-434,-17,0,0,0,68,-104,-6364,-12579,-5766,-117,0,0,0,0,0,-4,-57,-20,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-110,-1877,-356,-56,-10,0,0,-4,-44,-135,-1531,-8987,-8674,-951,-1,0,-2,3,-7,20,-7,-78,-4,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-778,-1602,-65,0,0,0,0,-520,-6437,-5314,-3547,-7959,-1541,-36,0,0,-274,88,-25,-181,-263,-30,0,0,0,-202,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,157,-2098,-235,-3,0,0,0,0,-6517,-21515,-3248,-3527,-3282,-67,0,0,0,-2891,-701,168,-772,-159,-2,0,0,0,-2135,-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,239,-456,-17,0,0,0,0,0,-9535,29066,62,-933,-175,-1,0,0,0,-3513,-2604,147,-221,-7,0,0,0,0,-2589,-62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,46,-1,0,0,0,0,0,-1622,10365,75,-35,-1,0,0,0,0,-491,-574,7,-6,0,0,0,0,0,-361,-8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-28,628,-1849,-745,0,0,0,0,0,-90,-542,-2782,-701,0,0,0,0,0,0,-10,-196,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-459,1134,-2482,-569,0,0,0,0,0,-1929,-6205,-3075,-1397,0,0,0,0,0,0,-58,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,-15,-1209,-1520,-386,-37,0,0,0,-13,-352,-8734,-12537,-867,-84,0,0,0,0,0,-7,-50,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,-335,-298,-619,-434,-17,0,0,0,68,-104,-6364,-12580,-5767,-117,0,0,0,0,0,-4,-57,-20,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-110,-1878,-356,-56,-10,0,0,-4,-44,-135,-1531,-8988,-8675,-951,-1,0,-2,3,-7,20,-7,-78,-4,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-778,-1602,-65,0,0,0,0,-520,-6439,-5314,-3547,-7960,-1541,-36,0,0,-274,88,-25,-181,-263,-30,0,0,0,-202,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,157,-2099,-235,-3,0,0,0,0,-6518,-21524,-3249,-3527,-3283,-67,0,0,0,-2891,-701,168,-772,-159,-2,0,0,0,-2135,-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,239,-456,-17,0,0,0,0,0,-9535,29062,62,-933,-175,-1,0,0,0,-3514,-2604,147,-221,-7,0,0,0,0,-2589,-62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,46,-1,0,0,0,0,0,-1622,10365,75,-35,-1,0,0,0,0,-492,-574,7,-6,0,0,0,0,0,-361,-8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-28,628,-1849,-745,0,0,0,0,0,-90,-542,-2782,-702,0,0,0,0,0,0,-10,-196,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-459,1134,-2482,-569,0,0,0,0,0,-1929,-6206,-3075,-1397,0,0,0,0,0,0,-58,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,-15,-1209,-1520,-386,-37,0,0,0,-13,-352,-8734,-12537,-867,-84,0,0,0,0,0,-7,-50,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,-335,-298,-619,-434,-17,0,0,0,68,-104,-6364,-12580,-5767,-117,0,0,0,0,0,-4,-57,-20,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-110,-1878,-356,-56,-10,0,0,-4,-44,-135,-1531,-8988,-8675,-951,-1,0,-2,3,-7,20,-7,-78,-4,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-778,-1602,-65,0,0,0,0,-520,-6439,-5315,-3547,-7961,-1541,-36,0,0,-274,88,-25,-181,-263,-30,0,0,0,-202,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,157,-2099,-235,-3,0,0,0,0,-6517,-21529,-3249,-3528,-3283,-67,0,0,0,-2892,-701,168,-772,-159,-2,0,0,0,-2135,-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,239,-456,-17,0,0,0,0,0,-9535,29057,62,-933,-175,-1,0,0,0,-3514,-2604,147,-221,-7,0,0,0,0,-2589,-62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,46,-1,0,0,0,0,0,-1622,10364,75,-35,-1,0,0,0,0,-492,-574,7,-6,0,0,0,0,0,-361,-8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-28,628,-1849,-745,0,0,0,0,0,-90,-542,-2782,-702,0,0,0,0,0,0,-10,-196,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-459,1134,-2482,-569,0,0,0,0,0,-1929,-6205,-3075,-1397,0,0,0,0,0,0,-58,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,-15,-1209,-1519,-386,-37,0,0,0,-13,-352,-8734,-12537,-867,-84,0,0,0,0,0,-7,-50,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,-335,-298,-619,-434,-17,0,0,0,68,-104,-6364,-12580,-5767,-117,0,0,0,0,0,-4,-57,-20,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-110,-1877,-356,-56,-10,0,0,-4,-44,-135,-1531,-8987,-8675,-951,-1,0,-2,3,-7,20,-7,-78,-4,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-778,-1602,-65,0,0,0,0,-520,-6439,-5314,-3547,-7960,-1541,-36,0,0,-274,88,-25,-181,-263,-30,0,0,0,-202,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,157,-2098,-235,-3,0,0,0,0,-6517,-21533,-3249,-3527,-3283,-67,0,0,0,-2891,-701,168,-772,-159,-2,0,0,0,-2135,-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5,239,-456,-17,0,0,0,0,0,-9534,29048,62,-933,-175,-1,0,0,0,-3513,-2604,147,-221,-7,0,0,0,0,-2589,-62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,46,-1,0,0,0,0,0,-1622,10362,75,-35,-1,0,0,0,0,-492,-574,7,-6,0,0,0,0,0,-361,-8,0,0,0,0,0,0,0};

//static float u0_max = 0, u1_max = 0, u0_min = 0, u1_min = 0;

void controllerGeomReset(void)
{
  t = 0;
}

void controllerGeomInit(void)
{
  controllerGeomReset();
/*  for(int i = 0; i < 2025; i++){
    if((float)ltb0[i]/100000.0f > u0_max) u0_max = (float)ltb0[i]/100000.0f;
    if((float)ltb1[i]/100000.0f > u1_max) u1_max = (float)ltb1[i]/100000.0f;
    if((float)ltb0[i]/100000.0f < u0_min) u0_min = (float)ltb0[i]/100000.0f;
    if((float)ltb1[i]/100000.0f < u0_min) u1_min = (float)ltb1[i]/100000.0f;
  }
*/}

bool controllerGeomTest(void)
{
  return true;
}

void controllerGeom(control_t *control, setpoint_t *setpoint,
                                         const sensorData_t *sensors,
                                         const state_t *state,
                                         const uint32_t tick)
{
  if (!RATE_DO_EXECUTE(ATTITUDE_RATE, tick)) {
      return;
    }

  struct vec er, ev, r1, r2, r3, M;
  struct mat33 Rd;
  float yaw_des = radians(setpoint->attitude.yaw);
  struct vec setpointPos = mkvec(setpoint->position.x, setpoint->position.y, setpoint->position.z);
  struct vec setpointVel = mkvec(setpoint->velocity.x, setpoint->velocity.y, setpoint->velocity.z);
  struct vec statePos = mkvec(state->position.x, state->position.y, state->position.z);
  struct vec stateVel = mkvec(state->velocity.x, state->velocity.y, state->velocity.z);
  float T;

  // Position Error (ep)
  er = vsub(statePos, setpointPos);

  // Velocity Error (ev)
  ev = vsub(stateVel, setpointVel);

  struct vec target_thrust = vzero();

  /* if (setpoint->mode.x == modeAbs) {
    target_thrust.x = -kr_xy*er.x - kv_xy*ev.x + g_vehicleMass*setpoint->acceleration.x;
    target_thrust.y = -kr_xy*er.y - kv_xy*ev.y + g_vehicleMass*setpoint->acceleration.y;
    target_thrust.z = -kr_z*er.z - kv_z*ev.z + g_vehicleMass*(setpoint->acceleration.x + GRAVITY_MAGNITUDE);
  } else {
    target_thrust.x = -sinf(radians(setpoint->attitude.pitch));
    target_thrust.y = -sinf(radians(setpoint->attitude.roll));
    if (setpoint->mode.z == modeAbs) {
      target_thrust.z = -kr_z*er.z - kv_z*ev.z + g_vehicleMass*GRAVITY_MAGNITUDE;
    } else {
      target_thrust.z = 1;
    }
  }*/
  
  if (!mode) {      // position control
    target_thrust.x = -kr_xy*er.x - kv_xy*ev.x + g_vehicleMass*setpoint->acceleration.x;
    target_thrust.y = -kr_xy*er.y - kv_xy*ev.y + g_vehicleMass*setpoint->acceleration.y;
    target_thrust.z = -kr_z*er.z - kv_z*ev.z + g_vehicleMass*(setpoint->acceleration.z + GRAVITY_MAGNITUDE);
    r3 = vnormalize(target_thrust);
    r2 = vnormalize(vcross(r3,mkvec(cosf(yaw_des), sinf(yaw_des), 0)));
    r1 = vcross(r2, r3);
    Rd = mcolumns(r1, r2, r3);
    wd = mkvec(radians(setpoint->attitudeRate.roll), -radians(setpoint->attitudeRate.pitch), radians(setpoint->attitudeRate.yaw));
  } else {
    target_thrust.z = -kr_z*er.z - kv_z*ev.z + g_vehicleMass*(setpoint->acceleration.z + GRAVITY_MAGNITUDE);
    t = getTime();
    T = t/timescale;
    // flip trajectory
    q0 = 1.998f/(1.0f+expf(-20.0f*(T-0.45f)))-0.999f;
    q2 = sqrtf(1.0f-q0*q0);
    float dq0 = -(39.9f * expf(-20.0f * (T - 0.45f))) / ( (expf(-20.0f * (T - 0.45f)) + 1.0f) * (expf(-20.0f * (T - 0.45f)) + 1.0f) );
    struct quat qd = mkquat(0, q2, 0, q0);
    struct vec eul_des = quat2rpy(qd);
    p_des = eul_des.y;
    Rd = quat2rotmat(qd);
    wd = mkvec(0.0f, -2.0f*dq0/sqrtf(1.0f-q0*q0), 0.0f);
  }
  q = mkquat(state->attitudeQuaternion.x, state->attitudeQuaternion.y, state->attitudeQuaternion.z, state->attitudeQuaternion.w);
  struct vec eul = quat2rpy(q);
  p = eul.y;
  struct mat33 R = quat2rotmat(q);

  struct mat33 eR1 = mmul(mtranspose(Rd),R);
  struct mat33 eR2 = mmul(mtranspose(R),Rd);
  struct mat33 tr_tmp = msub(meye(),eR1);
  psi = 0.5f*(tr_tmp.m[0][0] + tr_tmp.m[1][1] + tr_tmp.m[2][2]);

  struct mat33 eRm = msub(eR1,eR2);

  eR.x = 1.0f*eRm.m[2][1];
  eR.y = -1.0f*eRm.m[0][2];
  eR.z = 1.0f*eRm.m[1][0];

  float stateAttitudeRateRoll = radians(sensors->gyro.x);
  float stateAttitudeRatePitch = -radians(sensors->gyro.y);
  float stateAttitudeRateYaw = radians(sensors->gyro.z);

  struct vec w = mkvec(stateAttitudeRateRoll, stateAttitudeRatePitch, stateAttitudeRateYaw);

  // struct vec w_target = mvmul(eR2,wd);

  ew = vsub(w,wd);


  if (prev_omega_roll == prev_omega_roll) { /*d part initialized*/
    if (!mode){
      err_d_roll = ((radians(setpoint->attitudeRate.roll) - prev_setpoint_omega_roll) - (stateAttitudeRateRoll - prev_omega_roll)) / dt;
      err_d_pitch = (-(radians(setpoint->attitudeRate.pitch) - prev_setpoint_omega_pitch) - (stateAttitudeRatePitch - prev_omega_pitch)) / dt;
      prev_omega_roll = stateAttitudeRateRoll;
      prev_omega_pitch = stateAttitudeRatePitch;
      prev_setpoint_omega_roll = radians(setpoint->attitudeRate.roll);
      prev_setpoint_omega_pitch = radians(setpoint->attitudeRate.pitch);
    } else {
      err_d_roll = clamp(((wd.x - prev_setpoint_omega_roll) - (stateAttitudeRateRoll - prev_omega_roll)) / dt, -200, 200);
      err_d_pitch = clamp(((wd.y - prev_setpoint_omega_pitch) - (stateAttitudeRatePitch - prev_omega_pitch)) / dt, -200, 200);
      prev_omega_roll = stateAttitudeRateRoll;
      prev_omega_pitch = stateAttitudeRatePitch;
      prev_setpoint_omega_roll = wd.x;
      prev_setpoint_omega_pitch = wd.y;
    }
  }
  

  struct vec cross = vcross(w, mkvec(Ixx*w.x, Ixx*w.x, Izz*w.z));
  // cross = vzero();

  float thrust = 0;
  if (!mode){
    M.x = cross.x - kR_xy * eR.x - kw_xy * ew.x + kd_omega_rp * err_d_roll;
    M.y = cross.y - kR_xy * eR.y - kw_xy * ew.y + kd_omega_rp * err_d_pitch;
    M.z = -kR_z  * eR.z - kw_z  * ew.z;
    thrust = vdot(target_thrust, mcolumn(R,2));
  } else {
    /*float xtilde[4] = {q.x, q.y, w.x / 10.0f, w.y / 10.0f};
    float utmp = 0;

    float tmp1 = 0;
    float k1[48];
    for (int i = 0; i < 48; i++) {
      for (int j = 0; j < 4; j++) {
        tmp1 += -0.5f * lam1[j] * powf((xtilde[j] - X[i][j]), 2.0f);
      }
      k1[i] += sf1 * expf(tmp1);
    }
    for (int i = 0; i < 48; i++) {
      utmp += k1[i] * alpha1[i];
    }
    u1 = utmp / 200.0f;
    utmp = 0.0f;

    float tmp2 = 0;
    float k2[48];
    for (int i = 0; i < 48; i++) {
      for (int j = 0; j < 4; j++) {
        tmp2 += -0.5f * lam2[j] * powf((xtilde[j] - X[i][j]), 2.0f);
      }
      k2[i] += sf2 * expf(tmp2);
    }
    for (int i = 0; i < 48; i++) {
      utmp += k2[i] * alpha2[i];
    }
    u2 = utmp / 200.0f;
    */
    //if(gp){
    //  float xtilde[4] = {q.x, q.y, w.x / 10.0f, w.y / 10.0f};
      /*float min0 = 10;
      float min1 = 10;
      float min2 = 10;
      float min3 = 10;
      min0_i = 0; min1_i = 0; min2_i = 0; min3_i = 0;
      for(int i=0; i < 5; i++){
        float diff = fabsf(gr0[i] - xtilde[0]);
        if(diff < min0) {
          min0 = diff;
          min0_i = i;
        }
      }
      for(int i=0; i < 9; i++){
        float diff = fabsf(gr1[i] - xtilde[1]);
        if(diff < min1) {
          min1 = diff;
          min1_i = i;
        }
      }
      for(int i=0; i < 5; i++){
        float diff = fabsf(gr2[i] - xtilde[2]);
        if(diff < min2) {
          min2 = diff;
          min2_i = i;
        }
      }
      for(int i=0; i < 9; i++){
        float diff = fabsf(gr3[i] - xtilde[3]);
        if(diff < min3) {
          min3 = diff;
          min3_i = i;
        }
      }
      u1 = (float)ltb1[min0_i*min1_i*min2_i*min3_i] / 10000.0f / 200.0f / 10.0f;
      u0 = (float)ltb0[min0_i*min1_i*min2_i*min3_i] / 10000.0f / 200.0f / 10.0f;
      */
    //  u0 = LI_4D(5, 9, 5, 9, gr0, gr1, gr2, gr3, ltb0, xtilde[0], xtilde[1], xtilde[2], xtilde[3]);
    //  u0 = clamp(u0, u0_min, u0_max) / 200.0f;
    //  u1 = LI_4D(5, 9, 5, 9, gr0, gr1, gr2, gr3, ltb1, xtilde[0], xtilde[1], xtilde[2], xtilde[3]);
    //  u1 = clamp(u1, u1_min, u1_max) / 200.0f;
    //}
    
    M.x = cross.x - kR_xy * eR.x - kw_xy * ew.x + kd_omega_rp * err_d_roll;// - u0;
    M.y = cross.y - kR_xy * eR.y - kw_xy * ew.y + kd_omega_rp * err_d_pitch;// - u1;
    M.z = kR_z  * eR.z - kw_z  * ew.z;
    float den = R.m[2][2];
    if(den > 1e-7f){
      thrust_tmp = clamp(target_thrust.z / den, 0.22f, 0.5f);
    } else{
      thrust_tmp = 0.22f;
    }
    thrust = (0.22f * cosf(2.0f * 3.14f * t / (0.9f*timescale)) + 0.4f);
    //thrust = thrust_tmp;
  }
  
 
  if (setpoint->mode.z == modeDisable) {
    control->thrust = setpoint->thrust;
  } else {
    control->thrust = thrust * thrust_scale;
  }
  cmd_thrust_g = thrust*thrust_scale;
  r_roll = radians(sensors->gyro.x);
  r_pitch = -radians(sensors->gyro.y);
  r_yaw = radians(sensors->gyro.z);

  if (control->thrust > 0) {
    control->roll = clamp(M.x * thrust_scale / l, -32000, 32000);
    control->pitch = clamp(M.y * thrust_scale / l, -32000, 32000);
    control->yaw = clamp(-M.z * thrust_scale / b, -32000, 32000);
    cmd_roll = control->roll;
    cmd_pitch = control->pitch;
    cmd_yaw = control->yaw;
  } else {
    control->roll = 0;
    control->pitch = 0;
    control->yaw = 0;

    cmd_roll = control->roll;
    cmd_pitch = control->pitch;
    cmd_yaw = control->yaw;
  }
}

void setMode(bool val) {
  mode = val;
}

float getPsi(){
  return psi;
}

PARAM_GROUP_START(ctrlGeom)
PARAM_ADD(PARAM_FLOAT, kr_xy, &kr_xy)
PARAM_ADD(PARAM_FLOAT, kv_xy, &kv_xy)
PARAM_ADD(PARAM_FLOAT, kr_z, &kr_z)
PARAM_ADD(PARAM_FLOAT, kv_z, &kv_z)
PARAM_ADD(PARAM_FLOAT, kR_xy, &kR_xy)
PARAM_ADD(PARAM_FLOAT, kR_z, &kR_z)
PARAM_ADD(PARAM_FLOAT, kw_xy, &kw_xy)
PARAM_ADD(PARAM_FLOAT, kw_z, &kw_z)
PARAM_ADD(PARAM_FLOAT, kd_omega_rp, &kd_omega_rp)
PARAM_ADD(PARAM_FLOAT, timescale, &timescale)
PARAM_ADD(PARAM_UINT8, mode, &mode)
//PARAM_ADD(PARAM_UINT8, gp, &gp)
PARAM_ADD(PARAM_FLOAT, mass, &g_vehicleMass)
PARAM_GROUP_STOP(ctrlGeom)

LOG_GROUP_START(ctrlGeom)
LOG_ADD(LOG_FLOAT, cmd_thrust_g, &cmd_thrust_g)
LOG_ADD(LOG_FLOAT, cmd_roll, &cmd_roll)
LOG_ADD(LOG_FLOAT, cmd_pitch, &cmd_pitch)
LOG_ADD(LOG_FLOAT, cmd_yaw, &cmd_yaw)
LOG_ADD(LOG_FLOAT, r_roll, &r_roll)
LOG_ADD(LOG_FLOAT, r_pitch, &r_pitch)
//LOG_ADD(LOG_FLOAT, t, &t)
LOG_ADD(LOG_FLOAT, thrust, &thrust_tmp)
//LOG_ADD(LOG_FLOAT, p_des, &p_des)
//LOG_ADD(LOG_FLOAT, p, &p)
//LOG_ADD(LOG_FLOAT, wdy, &wd.y)
//LOG_ADD(LOG_FLOAT, err_d_pitch, &err_d_pitch)
LOG_ADD(LOG_FLOAT, qw, &q.w)
LOG_ADD(LOG_FLOAT, eRy, &eR.y)
LOG_ADD(LOG_FLOAT, ewy, &ew.y)
//LOG_ADD(LOG_FLOAT, u1, &u1)
//LOG_ADD(LOG_FLOAT, u0, &u0)
//LOG_ADD(LOG_FLOAT, psi, &psi)
// LOG_ADD(LOG_INT16, min0, &min0_i)
// LOG_ADD(LOG_INT16, min1, &min1_i)
// LOG_ADD(LOG_INT16, min2, &min2_i)
// LOG_ADD(LOG_INT16, min3, &min3_i)
// LOG_ADD(LOG_FLOAT, u2, &u2)
LOG_GROUP_STOP(ctrlGeom)


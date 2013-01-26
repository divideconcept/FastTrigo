#include "fasttrigo.h"

const float invtwopi=0.1591549f;
const float twopi=6.283185f;
const float threehalfpi=4.7123889f;
const float invhalfpi=0.6366198f;
const float pi=3.141593f;
const float halfpi=1.570796f;
const float quarterpi=0.7853982f;
static const __m128 SIGNMASK =  _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
//various ninja tricks:
//http://fastcpp.blogspot.fr/2011/03/changing-sign-of-float-values-using-sse.html
//http://www.songho.ca/misc/sse/sse.html
//http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/
//http://www.masmforum.com/board/index.php?PHPSESSID=786dd40408172108b65a5a36b09c88c0&topic=9515.0

//FastTrigo
namespace FastTrigo
{
    float atan(float x);
    float cos_s(float x);
};

//http://cbloomrants.blogspot.fr/2010/11/11-20-10-function-approximation-by_20.html
//http://assemblyrequired.crashworks.org/2009/10/16/timing-square-root/
//max error: 0.032% average error: 0.0094%
float FastTrigo::sqrt(float squared)
{
#ifdef FASTTRIGOACCURATE
    return _mm_cvtss_f32(_mm_sqrt_ss(_mm_set_ss(squared)));
#else
    return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(squared)))*squared;
#endif
}
float FastTrigo::length(float x, float y)
{
    return FastTrigo::sqrt(x*x+y*y);
}
float FastTrigo::length(float x, float y, float z)
{
    return FastTrigo::sqrt(x*x+y*y+z*z);
}


//http://nghiaho.com/?p=997
//http://www.researchgate.net/publication/3321724_Efficient_approximations_for_the_arctangent_function
//max error: 0.0015 radians (0.086 degrees) 0.024%
#ifndef FASTTRIGOACCURATE
float FastTrigo::atan(float x)
{
    return quarterpi*x - x*(fabs(x) - 1)*(0.2447f + 0.0663f*fabs(x));
}
#else
float FastTrigo::atan(float x)
{
    float u=x*x;
    float u2=u*u;
    float u3=u2*u;
    float u4=u3*u;
    float f=1.f+0.33288950512027f*u-0.08467922817644f*u2+0.03252232640125f*u3-0.00749305860992f*u4;
    return x/f;
}
#endif
float FastTrigo::atan2(float y, float x)
{
    if(fabs(x)>fabs(y)) {
        float atan=FastTrigo::atan(y/x);
        if(x>0.f)
            return atan;
        else
            return y>0.f?atan+pi:atan-pi;
    } else {
        float atan=FastTrigo::atan(x/y);
        if(x>0.f)
            return y>0.f?halfpi-atan:-halfpi-atan;
        else
            return y>0.f?halfpi+atan:-halfpi+atan;
    }
}

//http://www.ganssle.com/approx/approx.pdf
#ifndef FASTTRIGOACCURATE
//max error: 0.06%
float FastTrigo::cos_s(float x)
{
    const float c1= 0.99940307f;
    const float c2=-0.49558072f;
    const float c3= 0.03679168f;
    float x2;      // The input argument squared
    x2=x * x;
    return (c1 + x2*(c2 + c3 * x2));
}
#else
//max error: 0.0007%
float FastTrigo::cos_s(float x)
{
    const float c1= 0.9999932946f;
    const float c2=-0.4999124376f;
    const float c3= 0.0414877472f;
    const float c4=-0.0012712095f;
    float x2;      // The input argument squared
    x2=x*x;
    return (c1 + x2*(c2 + x2*(c3 + c4*x2)));
}
#endif
float FastTrigo::cos(float angle){
    //clamp to the range 0..2pi
    angle=angle-floorf(angle*invtwopi)*twopi;
    angle=angle>0.f?angle:-angle;

    if(angle<halfpi) return FastTrigo::cos_s(angle);
    if(angle<pi) return -FastTrigo::cos_s(pi-angle);
    if(angle<threehalfpi) return -FastTrigo::cos_s(angle-pi);
    return FastTrigo::cos_s(twopi-angle);
}
float FastTrigo::sin(float angle){
    return FastTrigo::cos(halfpi-angle);
}
void FastTrigo::sincos(float angle, float *sin, float *cos){
    //clamp to the range 0..2pi
    angle=angle-floorf(angle*invtwopi)*twopi;
    float sinmultiplier=angle>0.f&&angle<pi?1.f:-1.f;
    angle=angle>0.f?angle:-angle;

    if(angle<halfpi) {
        *cos=FastTrigo::cos_s(angle);
        *sin=sinmultiplier*FastTrigo::sqrt(1.f-*cos**cos);
        return;
    }
    if(angle<pi) {
        *cos=-FastTrigo::cos_s(pi-angle);
        *sin=sinmultiplier*FastTrigo::sqrt(1.f-*cos**cos);
        return;
    }
    if(angle<threehalfpi) {
        *cos=-FastTrigo::cos_s(angle-pi);
        *sin=sinmultiplier*FastTrigo::sqrt(1.f-*cos**cos);
        return;
    }
    *cos=FastTrigo::cos_s(twopi-angle);
    *sin=sinmultiplier*FastTrigo::sqrt(1.f-*cos**cos);
    return;
}

//FastTrigoQt
#ifdef QT_GUI_LIB
float FastTrigoQt::length(QVector2D vector){
    return FastTrigo::length(vector.x(),vector.y());
}
float FastTrigoQt::length(QVector3D vector){
    return FastTrigo::length(vector.x(),vector.y(),vector.z());
}
float FastTrigoQt::angle(QVector2D vector){
    return FastTrigo::atan2(vector.y(),vector.x());
}
float FastTrigoQt::azimuth(QVector3D vector){
    return FastTrigo::atan2(vector.y(),vector.x());
}
float FastTrigoQt::inclination(QVector3D vector){
    float projectedlength=FastTrigoQt::length(vector.toVector2D());
    return FastTrigo::atan2(vector.z(),projectedlength);
}
float FastTrigoQt::x(QVector2D vector){
    return vector.x()*FastTrigo::cos(vector.y());
}
float FastTrigoQt::y(QVector2D vector){
    return vector.x()*FastTrigo::sin(vector.y());
}
float FastTrigoQt::x(QVector3D vector){
    return vector.x()*FastTrigo::cos(vector.y())*FastTrigo::cos(vector.z());
}
float FastTrigoQt::y(QVector3D vector){
    return vector.x()*FastTrigo::sin(vector.y())*FastTrigo::cos(vector.z());
}
float FastTrigoQt::z(QVector3D vector){
    return vector.x()*FastTrigo::sin(vector.z());
}
QVector2D FastTrigoQt::cartesian2polar(QVector2D vector){
    return QVector2D(FastTrigoQt::length(vector),FastTrigoQt::angle(vector));
}
QVector2D FastTrigoQt::polar2cartesian(QVector2D vector){
    float sin,cos;
    FastTrigo::sincos(vector.y(),&sin,&cos);
    return QVector2D(vector.x()*cos,
                     vector.x()*sin);
}
QVector3D FastTrigoQt::cartesian2spherical(QVector3D vector){
    return QVector3D(FastTrigoQt::length(vector),FastTrigoQt::azimuth(vector),FastTrigoQt::inclination(vector));
}
QVector3D FastTrigoQt::spherical2cartesian(QVector3D vector){
    float azimuthsin, azimuthcos;
    FastTrigo::sincos(vector.y(),&azimuthsin,&azimuthcos);
    float inclinationsin, inclinationcos;
    FastTrigo::sincos(vector.z(),&inclinationsin,&inclinationcos);
    return QVector3D(vector.x()*azimuthcos*inclinationcos,
                     vector.x()*azimuthsin*inclinationcos,
                     vector.x()*inclinationsin);
}
#endif

//FastTrigoSSE
namespace FastTrigoSSE
{
    __m128 atan(__m128 x);
    __m128 cos_s(__m128 x);
};

__m128 FastTrigoSSE::sqrt(__m128 squared)
{
#ifdef FASTTRIGOACCURATE
    return _mm_sqrt_ps(squared);
#else
    return _mm_mul_ps(_mm_rsqrt_ps(squared),squared);
#endif
}
__m128 FastTrigoSSE::length(__m128 x, __m128 y)
{
    return FastTrigoSSE::sqrt(_mm_add_ps(_mm_mul_ps(x,x),_mm_mul_ps(y,y)));
}
__m128 FastTrigoSSE::length(__m128 x, __m128 y, __m128 z)
{
    return FastTrigoSSE::sqrt(_mm_add_ps(_mm_add_ps(_mm_mul_ps(x,x),_mm_mul_ps(y,y)),_mm_mul_ps(z,z)));
}

#ifndef FASTTRIGOACCURATE
__m128 FastTrigoSSE::atan(__m128 x)
{
    //                                      quarterpi*x
    //                                      - x*(fabs(x) - 1)
    //                                      *(0.2447f+0.0663f*fabs(x));
    return _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(quarterpi),x),
                      _mm_mul_ps(_mm_mul_ps(x,_mm_sub_ps(_mm_andnot_ps(SIGNMASK,x),_mm_set1_ps(1.f))),
                                 (_mm_add_ps(_mm_set1_ps(0.2447f),_mm_mul_ps(_mm_set1_ps(0.0663f),_mm_andnot_ps(SIGNMASK, x))))));
}
#else
__m128 FastTrigoSSE::atan(__m128 x)
{
    __m128 u=_mm_mul_ps(x,x);
    __m128 u2=_mm_mul_ps(u,u);
    __m128 u3=_mm_mul_ps(u2,u);
    __m128 u4=_mm_mul_ps(u3,u);
    //__m128 f=1.f+0.33288950512027f*u-0.08467922817644f*u2+0.03252232640125f*u3-0.00749305860992f*u4;

    __m128 f=_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_set1_ps(1.f),
    _mm_mul_ps(_mm_set1_ps(0.33288950512027f),u)),
    _mm_mul_ps(_mm_set1_ps(-0.08467922817644f),u2)),
    _mm_mul_ps(_mm_set1_ps(0.03252232640125f),u3)),
    _mm_mul_ps(_mm_set1_ps(-0.00749305860992f),u4));
    return _mm_div_ps(x,f);
}
#endif
__m128 FastTrigoSSE::atan2(__m128 y, __m128 x)
{
    __m128 absxgreaterthanabsy=_mm_cmpgt_ps(_mm_andnot_ps(SIGNMASK,x),_mm_andnot_ps(SIGNMASK,y));
    __m128 ratio=_mm_div_ps(_mm_add_ps(_mm_and_ps(absxgreaterthanabsy,y),_mm_andnot_ps(absxgreaterthanabsy,x)),
                            _mm_add_ps(_mm_and_ps(absxgreaterthanabsy,x),_mm_andnot_ps(absxgreaterthanabsy,y)));
    __m128 atan=FastTrigoSSE::atan(ratio);

    __m128 xgreaterthan0=_mm_cmpgt_ps(x,_mm_set1_ps(0.f));
    __m128 ygreaterthan0=_mm_cmpgt_ps(y,_mm_set1_ps(0.f));

    atan=_mm_xor_ps(atan,_mm_andnot_ps(absxgreaterthanabsy,_mm_and_ps(xgreaterthan0,SIGNMASK))); //negate atan if absx<=absy & x>0

    __m128 shift=_mm_set1_ps(pi);
    shift=_mm_sub_ps(shift,_mm_andnot_ps(absxgreaterthanabsy,_mm_set1_ps(halfpi))); //substract halfpi if absx<=absy
    shift=_mm_xor_ps(shift,_mm_andnot_ps(ygreaterthan0,SIGNMASK)); //negate shift if y<=0
    shift=_mm_andnot_ps(_mm_and_ps(absxgreaterthanabsy,xgreaterthan0),shift); //null if abs>absy & x>0

    return _mm_add_ps(atan,shift);
}

#ifndef FASTTRIGOACCURATE
__m128 FastTrigoSSE::cos_s(__m128 x)
{
    const __m128 c1=_mm_set1_ps( 0.99940307f);
    const __m128 c2=_mm_set1_ps(-0.49558072f);
    const __m128 c3=_mm_set1_ps( 0.03679168f);
    __m128 x2;      // The input argument squared
    x2=_mm_mul_ps(x,x);
    //               (c1+           x2*          (c2+           c3*x2));
    return _mm_add_ps(c1,_mm_mul_ps(x2,_mm_add_ps(c2,_mm_mul_ps(c3,x2))));
}
#else
__m128 FastTrigoSSE::cos_s(__m128 x)
{
    const __m128 c1=_mm_set1_ps( 0.9999932946f);
    const __m128 c2=_mm_set1_ps(-0.4999124376f);
    const __m128 c3=_mm_set1_ps( 0.0414877472f);
    const __m128 c4=_mm_set1_ps(-0.0012712095f);
    __m128 x2;      // The input argument squared
    x2=_mm_mul_ps(x,x);
    //               (c1+           x2*          (c2+           x2*          (c3+           c4*x2)));
    return _mm_add_ps(c1,_mm_mul_ps(x2,_mm_add_ps(c2,_mm_mul_ps(x2,_mm_add_ps(c3,_mm_mul_ps(c4,x2))))));
}
#endif
__m128 FastTrigoSSE::cos(__m128 angle){
    //clamp to the range 0..2pi

    //take absolute value
    angle=_mm_andnot_ps(SIGNMASK,angle);
    //fmod(angle,twopi)
    angle=_mm_sub_ps(angle,_mm_mul_ps(_mm_cvtepi32_ps(_mm_cvttps_epi32(_mm_mul_ps(angle,_mm_set1_ps(invtwopi)))),_mm_set1_ps(twopi))); //simplied SSE2 fmod, must always operate on absolute value
    //if SSE4.1 is always available, comment the line above and uncomment the line below
    //angle=_mm_sub_ps(angle,_mm_mul_ps(_mm_floor_ps(_mm_mul_ps(angle,_mm_set1_ps(invtwopi))),_mm_set1_ps(twopi))); //faster if SSE4.1 is always available

    __m128 cosangle=angle;
    cosangle=_mm_xor_ps(cosangle, _mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(halfpi)), _mm_xor_ps(cosangle,_mm_sub_ps(_mm_set1_ps(pi),angle))));
    cosangle=_mm_xor_ps(cosangle,_mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(pi)),SIGNMASK));
    cosangle=_mm_xor_ps(cosangle, _mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(threehalfpi)), _mm_xor_ps(cosangle,_mm_sub_ps(_mm_set1_ps(twopi),angle))));

    __m128 result=FastTrigoSSE::cos_s(cosangle);

    result=_mm_xor_ps(result,_mm_and_ps(_mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(halfpi)),_mm_cmplt_ps(angle,_mm_set1_ps(threehalfpi))),SIGNMASK));
    return result;
}
__m128 FastTrigoSSE::sin(__m128 angle){
    return FastTrigoSSE::cos(_mm_sub_ps(_mm_set1_ps(halfpi),angle));
}
void FastTrigoSSE::sincos(__m128 angle, __m128 *sin, __m128 *cos){
    //clamp to the range 0..2pi

    //fmod(angle,twopi)
    __m128i dividedangle=_mm_cvttps_epi32(_mm_mul_ps(angle,_mm_set1_ps(invtwopi)));
    dividedangle=_mm_sub_epi32(dividedangle,_mm_srli_epi32(dividedangle,31)); //substract -1 if negative
    angle=_mm_sub_ps(angle,_mm_mul_ps(_mm_cvtepi32_ps(dividedangle),_mm_set1_ps(twopi)));
    //if SSE4.1 is always available, comment the 3 lines above and uncomment the line below
    //angle=_mm_sub_ps(angle,_mm_mul_ps(_mm_floor_ps(_mm_mul_ps(angle,_mm_set1_ps(invtwopi))),_mm_set1_ps(twopi))); //faster if SSE4.1 is always available

    __m128 sinmultiplier=_mm_xor_ps(_mm_set1_ps(1.f),_mm_and_ps(_mm_cmplt_ps(angle,_mm_set1_ps(0.f)),_mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(-pi)),SIGNMASK)));

    //abs(angle)
    angle=_mm_andnot_ps(SIGNMASK,angle);

    __m128 cosangle=angle;
    cosangle=_mm_xor_ps(cosangle, _mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(halfpi)), _mm_xor_ps(cosangle,_mm_sub_ps(_mm_set1_ps(pi),angle))));
    cosangle=_mm_xor_ps(cosangle,_mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(pi)),SIGNMASK));
    cosangle=_mm_xor_ps(cosangle, _mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(threehalfpi)), _mm_xor_ps(cosangle,_mm_sub_ps(_mm_set1_ps(twopi),angle))));

    __m128 result=FastTrigoSSE::cos_s(cosangle);

    result=_mm_xor_ps(result,_mm_and_ps(_mm_and_ps(_mm_cmpge_ps(angle,_mm_set1_ps(halfpi)),_mm_cmplt_ps(angle,_mm_set1_ps(threehalfpi))),SIGNMASK));

    *sin=_mm_mul_ps(sinmultiplier,FastTrigoSSE::sqrt(_mm_sub_ps(_mm_set1_ps(1.f),_mm_mul_ps(result,result))));

    *cos=result;

    return;
}
void FastTrigoSSE::interleave(__m128 x0x1x2x3, __m128 y0y1y2y3, __m128 *x0y0x1y1, __m128 *x2y2x3y3)
{
    *x0y0x1y1=_mm_unpacklo_ps(x0x1x2x3,y0y1y2y3);
    *x2y2x3y3=_mm_unpackhi_ps(x0x1x2x3,y0y1y2y3);
}
void FastTrigoSSE::deinterleave(__m128 x0y0x1y1, __m128 x2y2x3y3, __m128 *x0x1x2x3, __m128 *y0y1y2y3)
{
    *x0x1x2x3=_mm_shuffle_ps(x0y0x1y1,x2y2x3y3,_MM_SHUFFLE(2,0,2,0));
    *y0y1y2y3=_mm_shuffle_ps(x0y0x1y1,x2y2x3y3,_MM_SHUFFLE(3,1,3,1));
}

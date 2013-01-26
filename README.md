FastTrigo 1.0 (c) 2013 Robin Lobel
=========

  Fast yet accurate trigonometric functions

  FastTrigo: best approximated algorithms found on the web (see fasttrigo.cpp for detailed references)
  FastTrigoQt: convenience functions if using QVector2D/QVector3D class from Qt, based on FastTrigo
  FastTrigoSSE: packed trigonometry ported from FastTrigo algorithms (using SSE, SSE2, SSE3)

  Default accuracy:
  
    FastTrigo::sqrt max error: 0.032% (average error: 0.0094%)
    FastTrigo::atan2 max error: 0.024% (0.0015 radians, 0.086 degrees)
    FastTrigo::cos max error: 0.06%
    FastTrigo::sin max error: 0.06%

  You can optionally define FASTTRIGOACCURATE to be even more accurate:

  result samples for sqrt:
  
    default:  4.27122 8.24475 0.58313  3.60468
    accurate: 4.272   8.24621 0.583095 3.60555
    standard: 4.272   8.24621 0.583095 3.60555
  result samples for atan2:
  
    default:  1.21309 -2.89625 2.60208 0.587799
    accurate: 1.21202 -2.89661 2.60118 0.588
    standard: 1.21203 -2.89661 2.60117 0.588003
  result samples for cos:
  
    default:  0.0706043 -0.144992 0.877807 -0.989482
    accurate: 0.0707372 -0.1455   0.877583 -0.989992
    standard: 0.0707372 -0.1455   0.877583 -0.989992

  Speed test with default accuracy (MSVC2012 x64):
  
    FastTrigo::sqrt speed up: x2.5 (from standard sqrt)
    FastTrigo::atan2 speed up: x2.3 (from standard atan2)
    FastTrigo::sin/cos speed up: x1.9 (from standard sin/cos)
    FastTrigo::sincos speed up: x2.3 (from standard sin+cos)
    FastTrigoSSE::sqrt speed up: x8 (from standard sqrt)
    FastTrigoSSE::atan2 speed up: x7.3 (from standard atan2)
    FastTrigoSSE::sin/cos speed up: x4.9 (from standard sin/cos)
    FastTrigoSSE::sincos speed up: x6.2 (from standard sin+cos)

  Speed test with FASTTRIGOACCURATE defined (MSVC2012 x64):
  
    FastTrigo::sqrt speed up: x1.5 (from standard sqrt)
    FastTrigo::atan2 speed up: x1.7 (from standard atan2)
    FastTrigo::sin/cos speed up: x1.6 (from standard sin/cos)
    FastTrigo::sincos speed up: x1.8 (from standard sin+cos)
    FastTrigoSSE::sqrt speed up: x4.9 (from standard sqrt)
    FastTrigoSSE::atan2 speed up: x5.2 (from standard atan2)
    FastTrigoSSE::sin/cos speed up: x4.3 (from standard sin/cos)
    FastTrigoSSE::sincos speed up: x5.2 (from standard sin+cos)

  Distributed under BSD License

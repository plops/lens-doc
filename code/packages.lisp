(defpackage :base
  (:use :cl)
  (:export 
   #:num #:vec #:with-arrays #:v #:copy-vec
   #:.+ #:.- #:.* #:./ #:dot #:norm #:normalize
   #:cross #:.s #:vx #:vy #:vz #:v-spherical
   #:check-unit-vector #:check-range
   #:req #:+pi+ #:+2*pi+ #:+pi/2+
   #:mat
   #:m
   #:rotation-matrix
   #:determinant
   #:transpose
   #:m*
   #:chop))

(defpackage :raster
  (:use :cl :base)
  (:export
   #:bfp-type
   #:raster-line
   #:raster-circle
   #:raster-disk
   #:raster-triangle))

(defpackage :raytrace
  (:use :cl :base)
  (:export))

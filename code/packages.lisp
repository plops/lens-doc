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
   #:chop
   #:*read-default-float-format*))

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
  (:export
   #:refract-plane
   #:intersect-plane
   #:intersect-sphere
   #:refract-objective-detection
   #:refract-objective-illumination
   #:refract-thin-lens
   #:aberrate-index-plane))

(defpackage :project
  (:use :cl :base :raytrace)
  (:export
   #:prepare-out-of-focus
   #:calc-periphery-point)) 

(defpackage :plot-macros
  (:use :cl :base)
  (:export
   #:with-asy
   #:asy
   #:coord
   #:line
   #:line-colored
   #:circle
   #:arc))
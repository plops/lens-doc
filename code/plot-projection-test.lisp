(require :raytrace)

(defpackage :plot-projection-test
  (:use :cl :base :raytrace :project
	:plot-macros))

(in-package :plot-projection-test)

(let* ((ne 1.33)
       (n 1.52)
       (na 1.38)
       (ftl 164.5)
       (mag 63s0)
       (f (/ 164.5 mag))
       (c (.s ne (v 0 0 0)))
       (s (.s ne (v -6 0 -7)))
       (r (* ne 1.2))
       (rbfp (* f na))
       (h .01)
       (d (+ (* ne h)
	     (* n (- f h)))))
  (format nil "~a~%" (list 'f f 'rbfp rbfp 'd d 'delta-d (- (* n f) d))))
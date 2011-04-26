#.(require :raytrace)

(defpackage :plot-projection-test
  (:use :cl :base :raytrace :project
	:plot-macros))

(in-package :plot-projection-test)

(let* ((ne 1.33)
       (n 1.52)
       (na 1.38)
       (ftl 164.5)
       (mag 63s0)
       (f (/ ftl mag))
       (c (.s ne (v 0 0 0)))
       (s (.s ne (v -6s-3 0 -7s-3)))
       (r (* ne 1.2s-3))
       (rbfp (* f na))
       (h .01)
       (d (+ (* ne h)
	     (* n (- f h))))
       (slip-z (* -1 n (- f h)))
       (slip-center (v 0 0 slip-z))
       (slip-normal (v 0 0 -1))
       (obj-center (v))
       (obj-normal (v 0 0 -1))
       (bfp-center (v 0 0 f))
       (bfp-normal obj-normal))
  (format t "~a~%" (list 'f f 'rbfp rbfp 'd d 'delta-d-um (* 1000 (- (* n f) d))))
  (with-asy "/dev/shm/projection-test.asy"
    (asy "import three;")
    (asy "size(100000,100000);")
    (asy "currentprojection=perspective(
camera=(0,15,0),up=(0,0,1),target=~a,
zoom=1,angle=0,autoadjust=false);" (coord (v 0 0 (- d))))
    #+nil     (asy "currentprojection=orthographic(
camera=(0,10000,0),up=(0,0,1),target=~a,showtarget=false,center=true);" 
	 (coord (v 0 0 (- d))))

    (line (v 0 0 -4) (v 0 0 10)) ; optical axis
    (line (v rbfp 0 f) (v 0 0 f)) ; bfp
    (asy "draw(circle(~a,~a));" (coord (v 0 0 f)) rbfp) ; bfp round
    (let* ((gauss-center (v 0 0 (* n f -1)))
	   (alpha (* (/ 180 +pi+) (asin (/ na n))))
	   (gauss-periphery (v 0 0 (- rbfp (* n f))))) ; transmissive part of gaussian sphere
      (asy "draw(arc(~a,~a,~a,0,~a,0),red+linewidth(300));"
	   (coord gauss-center) (* n f) (- alpha) alpha)
      (asy "draw(circle(~a,~a,(0,1,0)));" (coord gauss-center) (* n f)))
    (line (v 0 0 (- d)) (v 1 0 (- d))) ; aberrated focus
    (line (v 0 0 slip-z) (v 1 0 slip-z)) ; coverslip
    (let ((dz (v 0 0 (- d)))
	  (normal (v 0 1 0)))
      (asy "draw(circle(~a,~a,~a));" ;; out-of-focus
	   (coord (.+ s dz)) 
	   r 
	   (coord normal))
      (asy "draw(circle(~a,~a,~a));" ;; target
	   (coord (.+ c dz)) 
	   r 
	   (coord normal))
      (let ((c (.+ dz c))
	    (s (.+ dz s)))
	(multiple-value-bind (x1 y1 hx hy) ;; billboard
	    (prepare-out-of-focus c s :r r)
	  (let ((rays 9))
	    (dotimes (i rays)
	      (let* ((per (calc-periphery-point (* (/ +2*pi+ rays) i)
						;; point on out of focus nucleus
						x1 y1 hx hy s))
		     (hf (normalize (.- c per)))) ;; direction on cone
		(multiple-value-bind (dir! start! f!) ;; corrected positions
		    (aberrate-index-plane c hf slip-center slip-normal (/ ne n))
		  (multiple-value-bind (r e) 
		      (refract-objective-detection start! dir!
						   n f obj-center obj-normal)
		    #+nil(line start!
			       (.+ start! (.s 2s0 dir!)))
		    (let ((bfp (intersect-plane e r bfp-center bfp-normal))) 
		      (line per f! e bfp))))
		))))))))
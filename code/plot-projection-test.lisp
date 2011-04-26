#.(require :raytrace)

(defpackage :plot-projection-test
  (:use :cl :base :raytrace :project
	:plot-macros))

(in-package :plot-projection-test)

(let* ((ne 1.33)
       (n 1.52)
       (ftl 164.5)
       (mag 63.0)
       (na 1.38)
       (f (/ ftl mag))
       (c (.s ne (v 0 0 0)))
       (s (.s ne (v -5e-3 0 -7e-3)))
       (r (* ne 1.2e-3))
       (rbfp (* f na))
       (h .01)
       (d (+ (* ne h)
	     (* n (- f h))))
       (slip-z (* -1 n (- f h)))
       (slip-center (v 0 0 slip-z))
       (slip-normal (v 0 0 -1))
       (obj-center (v))
       (obj-normal (v 0 0 -1))
       (tube-center (v 0 0 (+ f ftl)))
       (tube-normal (v 0 0 1))
       (field .04)
       (field-cam (* mag field)))
  (format t "~a~%" (list 'f f 'rbfp rbfp 'd d 'delta-d-um (* 1000 (- (* n f) d))))
  (with-asy "/dev/shm/projection-test.asy"
    (asy "import three;")
    (asy "size(10000000,10000000);")
    #+nil (asy "currentprojection=perspective(
camera=(0,205,0),up=(0,0,1),target=~a,
zoom=1,angle=0,autoadjust=false);" (coord (v 0 0 (- d))))
     (asy "currentprojection=orthographic(
camera=(0,10000,0),up=(0,0,1),target=~a,showtarget=false,center=true);" 
	 (coord (v 0 0 (- d))))

    (line (v 0 0 -4) (v 0 0 (+ f ftl ftl))) ; optical axis
    (line (v rbfp 0 f) (v 0 0 f)) ; bfp
    (line (v rbfp 0 (+ f ftl)) (v 0 0 (+ f ftl)))		; tubelens
    (line (v field-cam 0 (+ f ftl ftl)) (v 0 0 (+ f ftl ftl)))		; camera
    (asy "draw(circle(~a,~a));" (coord (v 0 0 f)) rbfp) ; bfp round
    (let* ((gauss-center (v 0 0 (* n f -1)))
	   (alpha (* (/ 180 +pi+) (asin (/ na n))))) ; transmissive part of gaussian sphere
      (asy "draw(arc(~a,~a,~a,0,~a,0),red);"
	   (coord gauss-center) (* n f) (- alpha) alpha)
      (asy "draw(circle(~a,~a,(0,1,0)));" (coord gauss-center) (* n f)))
    (line (v 0 0 (- d)) (v field 0 (- d))) ; aberrated focus
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
		    (multiple-value-bind (rr ee) (refract-thin-lens e r ftl
								    tube-center tube-normal)
		      (line per f! e ee (.+ ee (.s (* 1.2 ftl) rr))))))))))))))
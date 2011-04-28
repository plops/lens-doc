#.(require :raytrace)

(defpackage :plot-projection-test
  (:use :cl :base :raytrace :project
	:plot-macros))

(in-package :plot-projection-test)

;; to export vector format use asy -render 0 <asy file>

;; out-of-focus target coverslip gauss-sphere bfp tubelens camera

;; bundle diameter 1.1mm, beam broadening half angle u=3.18 deg (for green)
;; microlenses 0.5mm x 0.5mm, f=23.5mm, broadening half angle 0.23 deg
;; rod 2.5mm x 2.5mm
;; magnified to 4mm x 4mm
;; MMA schlieren lens: 91.784 mm
;; relay lens 1 80.67 mm, relay lens 2 161.26mm
;; MEMI tubelens: 222.8 mm (zoom 1.221) .. 445.5 mm (zoom 2.441)



(defun plot-aberrated-microscope (&key (ne 1.33) (c (.s ne (v)))
					     (s (.s ne (v -5e-3 0 -7e-3)))
					     (r (* ne 1.2e-3)) (h .01))
  "Trace rays from periphery of out-of-focus nucleus at S with radius
R into back focal plane (MMA) and camera (LCoS). NE is the index of
the embedding medium (e.g. water). H is the depth of the embedding
medium. All dimensions in mm (if embedded multiply with ne)."
  (declare (type vec c s)
	   (type num ne r h))
  (let* ((n 1.52) ;; immersion index
	 (ftl 164.5) ;; focal length tubelens, 
	 (mag 63.0) ;; magnification
	 (na 1.38) ;; numerical aperture
	 (f (/ ftl mag))
	(rbfp (* f na))
	(d (+ (* ne h) ; position of aberrated focus
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
     ;(asy "size(10000000,10000000);")
     (asy "size(1000,1000);")
     #+nil (asy "currentprojection=perspective(
camera=(0,205,0),up=(0,0,1),target=~a,
zoom=1,angle=0,autoadjust=false);" (coord (v 0 0 (- d))))
     (asy "currentprojection=orthographic(
camera=(0,1000000,0),up=(0,0,1),target=~a,showtarget=true,center=false);" 
	  (coord (v)))
;     (asy "leftbutton=new string[] {\"shift\",\"rotate\",\"zoom\",\"pan\"};")
     (line (v 0 0 -4) (v 0 0 (+ f ftl ftl)))	  ; optical axis
     (line (v rbfp 0 f) (v 0 0 f))		  ; bfp
     (line (v rbfp 0 (+ f ftl)) (v 0 0 (+ f ftl))) ; tubelens
     (line (v field-cam 0 (+ f ftl ftl)) (v 0 0 (+ f ftl ftl)))	; camera
     (asy "draw(circle(~a,~a));" (coord (v 0 0 f)) rbfp) ; bfp round
     (let* ((gauss-center (v 0 0 (* n f -1)))
	    (alpha (* (/ 180 +pi+) (asin (/ na n))))) ; transmissive part of gaussian sphere
       (asy "draw(arc(~a,~a,~a,0,~a,0),red);"
	    (coord gauss-center) (* n f) (- alpha) alpha)
       (asy "draw(circle(~a,~a,(0,1,0)));" (coord gauss-center) (* n f)))
     (line (v 0 0 (- d)) (v field 0 (- d))) ; aberrated focus
     (line (v 0 0 slip-z) (v 1 0 slip-z))   ; coverslip
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
		       (line per f! e ee (.+ ee (.s (* 1.2 ftl) rr)))))))))))))))

#+nil
(project-through-aberrated-microscope)

;; same code as above, but no plot output. just an array of positions
;; 1) out-of-focus target coverslip behind-slip   0 1 2 3
;; 2) gauss-sphere bfp                            0 1 2 3 4 5
;; tubelens 
;; 3) infront-cam camera behind-cam               7 8 9
;; 4) everything but behind-cam                   0..9

(defun project-through-aberrated-microscope (&key (ne 1.33) (c (.s ne (v)))
					     (s (.s ne (v -5e-3 0 -7e-3)))
					     (r (* ne 1.2e-3)) (h .01) (rays 9))
  "Trace rays from periphery of out-of-focus nucleus at S with radius
R into back focal plane (MMA) and camera (LCoS). NE is the index of
the embedding medium (e.g. water). H is the depth of the embedding
medium. All dimensions in mm (if embedded multiply with ne)."
  (declare (type vec c s)
	   (type num ne r h))
  (let* ((n 1.52) ;; immersion index
	 (ftl 164.5) ;; focal length tubelens
	 (mag 63.0) ;; magnification
	 (na 1.38) ;; numerical aperture
	 (f (/ ftl mag))
	 (rbfp (* f na))
	 (d (+ (* ne h) ; position of aberrated focus
	       (* n (- f h))))
	 (slip-z (* -1 n (- f h)))
	 (slip-center (v 0 0 slip-z))
	 (slip-behind (v 0 0 (* -1 n (- f .05)))) ;; for plotting part of the ray in glass
	 (slip-normal (v 0 0 -1))
	 (obj-center (v))
	 (obj-normal (v 0 0 -1))
	 (bfp-center (v 0 0 f))
	 (bfp-normal obj-normal)
	 (tube-center (v 0 0 (+ f ftl)))
	 (tube-normal (v 0 0 1))
	 (dz-cam 1e0)
	 (camera-infront (v 0 0 (+ f ftl ftl (- dz-cam)))) ;; for plotting rays around cam
	 (camera-center (v 0 0 (+ f ftl ftl)))
	 (camera-behind (v 0 0 (+ f ftl ftl (+ dz-cam))))  ;; for plotting rays around cam
	 (camera-normal obj-normal)
	 (field .04)
	 (field-cam (* mag field))
	 (result (make-array (list rays 10) :element-type 'vec))
	 (dz (v 0 0 (- d)))
	 (c (.+ dz c))
	 (s (.+ dz s)))
    (multiple-value-bind (x1 y1 hx hy) ;; billboard
	(prepare-out-of-focus c s :r r)
      (dotimes (i rays)
	(let* ((per (calc-periphery-point (* (/ +2*pi+ rays) i)
					  ;; point on out-of-focus nucleus
					  x1 y1 hx hy s))
	       (hf (normalize (.- c per)))) ;; direction on cone
	  (setf (aref result i 0) per) ;; point on out-of-focus nucleus
	  (setf (aref result i 1) c) ;; target
	  (multiple-value-bind (dir! start! f!) ;; corrected positions
	      (aberrate-index-plane c hf slip-center slip-normal (/ ne n))
	    (setf (aref result i 2) f!  ;; intersection with coverslip
		  (aref result i 3) ;; intersection somewhere in immersion, for plotting 
		  (intersect-plane start! dir! slip-behind slip-normal))
	    (multiple-value-bind (r e) 
		(refract-objective-detection start! dir!
					     n f obj-center obj-normal)
	      (setf (aref result i 4) e ;; intersection with gauss-sphere
		    (aref result i 5) ;; intersection with bfp 
		    (intersect-plane e r bfp-center bfp-normal))
	      (multiple-value-bind (rr ee) (refract-thin-lens e r ftl
							      tube-center tube-normal)
		(setf (aref result i 6) ee ;; intersection with tubelens
		      (aref result i 7) ;; for plotting
		      (intersect-plane ee rr camera-infront camera-normal)
		      (aref result i 8) ;; intersection with camera
		      (intersect-plane ee rr camera-center camera-normal)
		      (aref result i 9) ;; for plotting
		      (intersect-plane ee rr camera-behind camera-normal)))))))
      result)))
#+nil
(format t "~a~%"  (project-through-aberrated-microscope :rays 1))


(defun plot-rays-around-focus (fn ray-arrays)
  (declare (type (simple-array vec 2) ray-arrays)
	   (type string fn))
  (with-asy fn
    (asy "import three;")
    (asy "size(1000,1000);")
    (asy "currentprojection=orthographic(
camera=(0,1000000,0),up=(0,0,1),target=~a,showtarget=true,center=false);" 
	 (coord (aref ray-arrays 0 1)))
    (let ((rays (array-dimension ray-arrays 0)))
      (dotimes (i rays)
	(macrolet ((r (point-index)
		      `(aref ray-arrays i ,point-index)))
	  (line (r 0) (r 1) (r 2) (r 3)))))))

#+nil 
(plot-rays-around-focus "/dev/shm/focus.asy"
 (project-through-aberrated-microscope :rays 13))

(defun plot-rays-around-camera (fn ray-arrays)
  (declare (type (simple-array vec 2) ray-arrays)
	   (type string fn))
  (with-asy fn
    (asy "import three;")
    (asy "size(1000,1000);")
    (asy "currentprojection=orthographic(
camera=(0,1000000,0),up=(0,0,1),target=~a,showtarget=true,center=false);" 
	 (coord (aref ray-arrays 0 8)))
    (let ((rays (array-dimension ray-arrays 0)))
      (dotimes (i rays)
	(macrolet ((r (point-index)
		      `(aref ray-arrays i ,point-index)))
	  (line (r 7) (r 8) (r 9)))))))

#+nil
(plot-rays-around-camera "/dev/shm/camera.asy"
 (project-through-aberrated-microscope :rays 13))

(defun plot-rays-to-bfp (fn ray-arrays)
  (declare (type (simple-array vec 2) ray-arrays)
	   (type string fn))
  (with-asy fn
    (asy "import three;")
    (asy "size(1000,1000);")
    (asy "currentprojection=orthographic(
camera=(0,10,0),up=(0,0,1),target=~a,showtarget=true,center=false);" 
	 (coord (aref ray-arrays 0 4)))
    (let ((rays (array-dimension ray-arrays 0)))
      (dotimes (i rays)
	(macrolet ((r (point-index)
		      `(aref ray-arrays i ,point-index)))
	  (line (r 3) (r 4) (r 5)))))))

#+nil
(plot-rays-to-bfp "/dev/shm/bfp.asy"
 (project-through-aberrated-microscope :rays 13))
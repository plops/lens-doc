(require :raytrace)

(defpackage :project
  (:use :cl :base :raytrace))

(in-package :project)

(declaim (optimize (debug 3) (speed 1) (safety 3)))

(defmacro with-asy (fn &body body)
  (let ((s (gensym)))
    `(with-open-file (,s ,fn
			 :direction :output
			 :if-exists :supersede
			 :if-does-not-exist :create)
       (macrolet ((asy (str &rest rest)
		    `(progn
		       (format ,',s ,str ,@rest)
		       (terpri ,',s))))
	 (labels ((coord (v)
		    (format nil "(~f,~f,~f)" (vx v) (vy v) (vz v)))
		  (line (&rest vecs)
		    (when (listp vecs)
		      (format ,s "draw(~a" (coord (first vecs)))
		      (dolist (v (cdr vecs))
			(format ,s "--~a" (coord v)))
		      (format ,s ");~%"))))
	   ,@body)))))



(defun prepare-out-of-focus (target s &key (ne 1.33s0) (r 1.2s0))
  "Input dimensions are in um and not corrected for embedding index."
  (declare (type vec target s)
	   (type num r ne))
  (let* ((s (.s ne s))
	 (r (* ne r))
	 (c (.s ne target))
	 (x (.- c s))
	 (rr (norm x))
	 (q (if (< (vz x) (* 2/3 rr))
		(v 0 0 1)
		(v 0 1 0)))
	 (y (cross x q))
	 (hx (if (< rr 1e-5)
		 (error "the nucleus is too close")
		 (normalize x)))
	 (hy (normalize y))
	 (x1 (* r r 1/2 (/ rr)))
	 (y1 (if (< (* 2 rr) r)
		 (error "the nucleus is too close")
		 (* r 1/2 (/ rr) (sqrt (- (* 4 rr rr) (* r r)))))))
    (values x1 y1 hx hy s c r)))


(defun calc-e (phi x1 y1 hx hy s)
  "Sample the periphery of the bill board circle that is formed by the
touching cone on a nucleus."
  (declare (type num x1 y1 phi)
	   (type vec hx hy s))
  (the vec
   (.+ s 
       (.s x1 hx)
       (.s y1 (m* (rotation-matrix phi hx) hy)))))

(with-asy "/dev/shm/o.asy"
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  
  (let ((rays 90)
	(ne 1.3s0))
    (multiple-value-bind (x1 y1 hx hy s c r) (prepare-out-of-focus (v 0 0 .1) (v 2 0 2) :ne ne) 
      ;; green nucleus
      (asy "draw(shift(~a)*scale3(~f)*unitsphere,green);" (coord s) r)
      ;; lines on surface of cone
      (dotimes (i rays)
	(let* ((e (calc-e (* (/ +2*PI+ rays) i) x1 y1 hx hy s))
	       (f (.- c e))
	       (hf (normalize f))
	       (normal (v 0 0 1))
	       (i (intersect-plane e hf (v) normal))
	       (ht (refract-plane hf normal (/ ne 1.5))))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord e)
	       (coord i)
	       (coord (.+ i (.s 5s0 ht)))))))))

;; normal is directed towards sample (opposite i)
(defun refract-objective-detection (p i n f center normal)
  (declare (type vec p i normal center)
	   (type num n f))
  (let* ((gauss-center (.+ center (.s (* n f) normal)))
	(nf-hit (intersect-plane p i gauss-center  normal))
	(e (intersect-sphere p i gauss-center (* n f)))
	(a (.s (* f (- 1 n)) normal))
	(f-mid (.- nf-hit a))
	(m (normalize (.- center f-mid))))
    (values m e)))

#+nil
(refract-objective-detection (v .1 0 -3.1) (v 0 0 1) 1.5s0 2s0 (v) (v 0 0 -1))


(with-asy "/dev/shm/objective.asy"
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let ((rays 13)
	(p (v .1 0 -3.1)))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (r e) (refract-objective-detection p dir
								1.5 2.0 (v) (v 0 0 -1))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord p)
	       (coord e)
	       (coord (.+ e (.s 1s0 r)))))))))

;; normal is directed in the opposite direction direction of i
(defun refract-thin-lens (start i f center normal)
  (declare (type vec start i normal center)
	   (type num f))
  (let* ((lens-hit (intersect-plane start i center normal))
	 (rho (.- lens-hit center))
	 (cos-theta (dot i normal))
	 (r (.- (.s (/ f cos-theta) i)
		rho)))
    (values (normalize r) lens-hit)))

#+nil
(refract-thin-lens (v -.1 0 -2) (normalize (v .1 0 1))
		   2.0 (v) (v 0 0 1))

(with-asy "/dev/shm/tubelens.asy"
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let ((rays 13)
	(p (v -.1 0 -2)))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (r e) (refract-thin-lens p dir 2.0 (v) (v 0 0 1))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord p)
	       (coord e)
	       (coord (.+ e (.s 1s0 r)))))))))


(with-asy "/dev/shm/microscope.asy"
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let* ((rays 4)
	 (p (v .1 0 -3.1))
	 (obj-c (v))
	 (obj-n (v 0 0 -1))
	 (obj-f 2.0)
	 (tube-f 10.0)
	 (tube-c (v 0 0 (+ obj-f tube-f)))
	 (tube-n (v 0 0 1)))
    (dotimes (i rays)
      (dotimes (j rays)
       (let* ((dir (normalize (v (/ (- i (floor rays 2)) rays)
				 (/ (- j (floor rays 2)) rays)
				 1))))
	 (multiple-value-bind (r e) (refract-objective-detection p dir
								 1.5 obj-f obj-c obj-n)
	   (multiple-value-bind (rr ee) (refract-thin-lens e r tube-f tube-c tube-n)
	     (asy "draw((~a)--(~a)--(~a)--(~a));" 
		  (coord p)
		  (coord e)
		  (coord ee)
		  (coord (.+ ee (.s (* 1.2 tube-f) rr)))
		  ))))))))

;; normal opposing dir
(defun aberrate-index-plane (start dir center normal ne/n)
  "A light source is embedded in index ne. CENTER and NORMAL define a
plane interface to a region of index n. Find equivalent source
position and direction if index jump wouldn't exist."
  (declare (type vec start dir center normal)
	   (type num ne/n))
  (let* ((f (intersect-plane start dir center normal))
	 (dir! (refract-plane dir normal ne/n))
	 (scaled-photon-distance (* ne/n (norm (.- f start))))
	 (start! (.- f (.s scaled-photon-distance dir!))))
    (values dir! start! f)))

(with-asy "/dev/shm/objective-aberrate.asy"
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let ((rays 13)
	(p (v .1 0 -3.1))
	(slide-center (v 0 0 -1))
	(slide-normal (v 0 0 -1))
	(n 1.5)
	(ne 1.3))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (dir! start! f!)
	    (aberrate-index-plane p dir slide-center slide-normal (/ ne n))
	  (multiple-value-bind (r e) 
	      (refract-objective-detection start! dir!
					   n 2.0 (v) (v 0 0 -1))
	    (line p
		  f!
		  e
		  (.+ e (.s 1s0 r)))))))))
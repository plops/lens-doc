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
	 (flet ((coord (v)
		  (format nil "(~f,~f,~f)" (vx v) (vy v) (vz v))))
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
  
  (let ((rays 13))
    (multiple-value-bind (x1 y1 hx hy s c r) (prepare-out-of-focus (v 0 0 .1) (v 2 0 2)) 
      ;; green nucleus
      (asy "draw(shift(~a)*scale3(~f)*unitsphere,green);" (coord s) r)
      ;; lines on surface of cone
      (dotimes (i rays)
	(let* ((e (calc-e (* (/ +2*PI+ rays) i) x1 y1 hx hy s))
	       (f (.- c e))
	       (hf (normalize f))
	       (i (intersect-plane e hf (v) (v 0 0 1))))
	  (asy "draw((~a)--(~a));" 
	       (coord e)
	       (coord i)))))))
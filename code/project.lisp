(require :raytrace)

(defpackage :project
  (:use :cl :base :raytrace))

(in-package :project)

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

(let* ((ne 1.33s0)
       (c (v))
       (s (.s ne (v 6 0 7)))
       (r (* ne 1.2))
       (x (.- c s))
       (rr (norm x))
       (q (if (< (vz x) (* 2/3 rr))
	      (v 0 0 1)
	      (v 0 1 0)))
       (y (cross x q))
       (hx (normalize x))
       (hy (normalize y))
       (x1 (* r r 1/2 (/ rr)))
       (y1 (if (< rr r)
	       (error "the nucleus is too close")
	       (* r 1/2 (/ rr) (sqrt (- (* 4 rr rr) (* r r))))))
       (e (.+ s 
	      (.s x1 hx)
	      (.s y1 hy))))
 (with-asy "/dev/shm/o.asy"
   (asy "import three;")
   (asy "size(1000,1000);")
   ;; coordinate axes
   (asy "draw((0,0,0)--(1,0,0),red);")
   (asy "draw((0,0,0)--(0,1,0),green);")
   (asy "draw((0,0,0)--(0,0,1),blue);")
   ;; green nucleus
   (asy "draw(shift(~a)*scale3(~f)*unitsphere,green);" (coord s) r)
   ;; line of cone
   (asy "draw((~a)--(~a));" (coord e) (coord (v)))))
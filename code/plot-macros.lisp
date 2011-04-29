(in-package :plot-macros)

(defmacro with-asy (fn &body body)
  "Open the file FN and write asymptote commands into it. This is
useful for 3D drawings."
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
		      (format ,s ");~%")))
		  (line-colored (color &rest vecs) ;; FIXME maybe introduce draw instead
		    (when (listp vecs)
		      (format ,s "draw(~a" (coord (first vecs)))
		      (dolist (v (cdr vecs))
			(format ,s "--~a" (coord v)))
		      (format ,s ",~a);~%" color)))
		  (circle (center &optional (radius 1e0) (normal (v 0 0 1)))
		    (asy "draw(circle(~a,~a,~a));"
			 (coord center) radius (coord normal)))
		  (arc (center radius &key (theta1 -1e0) (phi1 0e0) (theta2 (- theta1)) (phi2 0e0) (normal (v 0 0 1)))
		    (let ((s (/ 180e0 +pi+)))
		     (asy "draw(arc(~a,~a,~a,~a,~a,~a,~a));"
			  (coord center) radius (* s theta1) (* s phi1) (* s theta2) (* s phi2) (coord normal)))))
	   ,@body)))))
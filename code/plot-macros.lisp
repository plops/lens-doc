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
		  (line-colored (color &rest vecs)
		    (when (listp vecs)
		      (format ,s "draw(~a" (coord (first vecs)))
		      (dolist (v (cdr vecs))
			(format ,s "--~a" (coord v)))
		      (format ,s ",~a);~%" color))))
	   ,@body)))))
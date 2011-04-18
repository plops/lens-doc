(in-package :raytrace)

(defun refract-plane (incident normal n1/n2)
  (declare (type vec incident normal)
	   (type num n1/n2))
  (let* ((in (dot incident normal))
	 (in2 (expt in 2))
	 (eta2 (expt n1/n2 2))
	 (rad (- 1 (* eta2 (- 1 in2)))))
    (the vec
     (if (< rad 0)
	 (.- incident
	     (.s (* 2 in) normal))
	 (.- (.s n1/n2 incident)
	     (.s (+ (* n1/n2 in) (sqrt rad)) normal))))))


#+nil
(refract-plane (v-spherical .1s0 0s0) (v 0 0 -1) (/ 1.5))


(defun intersect-plane (start dir c n)
  "Intersection of ray and plane"
  (declare (type vec start dir c n))
  (let* ((hesse-dist (dot c n))
         (div (dot n dir)))
    (when (< (abs div) 1s-5)
      (error "ray and plane are parallel"))
    (let ((tau (/ (- hesse-dist (dot n start))
                div)))
      (the vec (.+ start (.s tau dir))))))
#+nil
(intersect-plane (v -1 1) (v 1) (v) (v 1))




(defun quadratic-root (a b c)
  (declare (type num a b c)
	   (values num num &optional))
  (let ((eps 1s-7)
	(det2 (- (* b b) (* 4 a c))))
    (the num 
      (if (<= det2 0)
	  (error "no solution")
	  (let* ((q (* -.5s0 (+ b (* (signum b) (sqrt det2)))))
		 (aa (abs a))
		 (aq (abs q)))
	    (cond ((< aq eps) (values (/ q a) (/ q a)))
		  ((< aa eps) (values (/ c q) (/ c q)))
		  (t (values (/ q a) (/ c q)))))))))

#+nil
(let ((w '((0 1 -1 1)
	   (0 2 -1 .5)
	   (1 0 -4 (2 2))
	   (1 2 -3 (-3 1))
	   (3 4 -12 (-2.77 1.44)))))
  (dolist (e w)
    (destructuring-bind (a b c r) e
	(format t "~a~%"
		(list r 
		      (multiple-value-list
		       (quadratic-root (* 1s0 a) (* 1s0 b) (* 1s0 c))))))))

(defun intersect-sphere (start dir center radius)
  (declare (type vec start dir center)
           (type num radius))
 ; (check-unit-vector dir)
  (let* ((l (.- start center))
         (c (- (dot l l) (* radius radius)))
         (b (* 2s0 (dot l dir))))
    (the vec (multiple-value-bind (t1 t2)
		 (quadratic-root 1s0 b c)
	       ;; we are only interested in the intersection in
	       ;; forward direction of the ray. if there are two
	       ;; intersections in forward direction, return the
	       ;; closest
	       (.+ start (.s (cond ((< 0 (* t1 t2)) (min t1 t2))
				   (t (max t1 t2)))
			     dir))))))

#+nil
(intersect-sphere (v .1) (normalize (v 0 0 1)) (v) 2s0)
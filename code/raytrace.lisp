(in-package :raytrace)

(defun refract-plane (incident normal n1/n2)
  ;; direct normal opposite of incident (like for a mirror)
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
(refract-plane (v-spherical .1e0 0e0) (v 0 0 -1) (/ 1.5))



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
	  (let* ((q (* -.5e0 (+ b (* (signum b) (sqrt det2)))))
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
		       (quadratic-root (* 1e0 a) (* 1e0 b) (* 1e0 c))))))))



(defun intersect-sphere (start dir center radius)
  (declare (type vec start dir center)
           (type num radius))
 ; (check-unit-vector dir)
  (let* ((l (.- start center))
         (c (- (dot l l) (* radius radius)))
         (b (* 2e0 (dot l dir))))
    (the vec (multiple-value-bind (t1 t2)
		 (quadratic-root 1e0 b c)
	       ;; we are only interested in the intersection in
	       ;; forward direction of the ray. if there are two
	       ;; intersections in forward direction, return the
	       ;; closest
	       (.+ start (.s (cond ((< 0 (* t1 t2)) (min t1 t2))
				   (t (max t1 t2)))
			     dir))))))

#+nil
(intersect-sphere (v .1) (normalize (v 0 0 1)) (v) 2e0)




(defun refract-objective-detection (p i n f center normal)
  ;; normal is directed towards sample (opposite i)
  "Propagate ray defined by start point p and direction i from sample
through an objective defined by its focal length f, immersion index n
the center position CENTER (the position where the gaussian sphere
intersects the optical axis) and orientation given by NORMAL. NORMAL
is directed nearly opposite to i."
  (declare (type vec p i normal center) 
	   (type num n f))
  (let* ((gauss-center (.+ center (.s (* n f) normal)))
	 (nf-hit (intersect-plane p i gauss-center normal))
	 (e (intersect-sphere p i gauss-center (* n f)))
	 (a (.s (+ (* f (- n 1))) normal))
	 (f-mid (.- nf-hit a))
	 (m (normalize (.- center f-mid))))
    (values m e)))

#+nil
(refract-objective-detection (v .1 0 -3.1) (v 0 0 1) 1.5e0 2e0 (v) (v 0 0 -1))



(defun refract-objective-illumination (start dir n f center normal)
  ;; normal is directed in the direction of i
  "Propagate ray defined by start point START and direction DIR from
behind the objective into the sample. The objective is defined by its
focal length f, immersion index n the center position CENTER and
orientation given by NORMAL. NORMAL is directed in nearly the same
direction as DIR."
  (declare (type vec start dir normal center)
	   (type num f n))
  (let* ((lens-hit (intersect-plane start dir center normal))
	 (rho (.- lens-hit center))
	 (cos-theta (dot dir normal))
	 (r (.- (.s (/ f cos-theta) dir)
		rho))
	 (a (.s (* f (- n 1)) normal))
	 (nf (* n f))
	 ;; (rho2 (dot rho rho))
	 ;; (nf2 (* nf nf))
	 ;; (s (.s (- nf (sqrt (- nf2 rho2) )) dir)) ;; approximation
	 (gauss-center (.+ center (.s nf normal)))
	 (s-global (intersect-sphere start dir gauss-center nf))
	 (s (.- s-global lens-hit)))
    (values (.+ r (.- a s)) (.+ s lens-hit))))



(defun refract-thin-lens (start i f center normal)
  "Refract a ray defined by postion START and direction i through a
lens of focal length f. The lens position is defined by CENTER and its
orientation by NORMAL. START can be on either side of the lens."
  ;; normal is directed in the opposite direction of i
  (declare (type vec start i normal center)
	   (type num f))
  (let* ((lens-hit (intersect-plane start i center normal))
	 (rho (.- lens-hit center))
	 (cos-theta (dot i normal)) ;; FIXME does it matter if normal is inverted?
	 (r (.- (.s (/ f cos-theta) i)
		rho)))
    (values (normalize r) lens-hit)))

#+nil
(refract-thin-lens (v -.1 0 -2) (normalize (v .1 0 1))
		   2.0 (v) (v 0 0 1))



(defun aberrate-index-plane (start dir center normal ne/n)
  ;; normal opposing dir
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

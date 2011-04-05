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

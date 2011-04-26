(in-package :project)

(declaim (optimize (debug 3) (speed 1) (safety 3)))

(defun prepare-out-of-focus (target s &key (r 1.2s-3))
  "Input dimensions are in mm and corrected for embedding index."
  (declare (type vec target s)
	   (type num r)
	   (values num num vec vec &optional))
  (let* ((x (.- target s))
	 (rr (norm x))
	 (q (if (< (vz x) (* 2/3 rr))
		(v 0 0 1)
		(v 0 1 0)))
	 (y (cross x q))
	 (hx (if (< rr 1e-8)
		 (error "the nucleus is too close")
		 (normalize x)))
	 (hy (normalize y))
	 (x1 (* r r 1/2 (/ rr)))
	 (y1 (if (< (* 2 rr) r)
		 (error "the nucleus is too close")
		 (* r 1/2 (/ rr) (sqrt (- (* 4 rr rr) (* r r)))))))
    (values x1 y1 hx hy)))


(defun calc-periphery-point (phi x1 y1 hx hy s)
  "Sample the periphery of the bill board circle that is formed by the
touching cone on a nucleus."
  (declare (type num x1 y1 phi)
	   (type vec hx hy s))
  (the vec
   (.+ s 
       (.s x1 hx)
       (.s y1 (m* (rotation-matrix phi hx) hy)))))











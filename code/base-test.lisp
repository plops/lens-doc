(defpackage :base-test
  (:use :cl :base))
(in-package :base-test)

#+nil
(define-test cross
  (assert-equal (v 0 0 1) (cross (v 1) (v 0 1))))

(defun num~ (a b &optional (abserr 1e-2) (relerr 1e-2))
  "Compare two numbers, check that difference is smaller than absolute
error or relative error"
  (declare (type num a b))
  (let ((aa (abs a)))
    (when (or (< aa abserr)
	      (< (abs (- aa (abs b))) (* relerr aa)))
      t)))

(defun v= (a b)
  (declare (type vec a b))
  #+nil (the boolean
	  (< (norm (.- a b)) 1s-3))
  (dotimes (i 3)
    (unless (num~ (aref a i) (aref b i))
      (error "vectors not equal"))))

(defun vr (&optional (x 10s0))
  (v (- (random (* 2 x)) x)
     (- (random (* 2 x)) x)
     (- (random (* 2 x)) x)))

;; some cross product identities from mathworld

(defun c10 (a b)
  (declare (type vec a b))
  (v= (cross a b)
      (.s -1s0 (cross b a))))

(defun c11 (a b c)
  (declare (type vec a b c))
  (v= (cross a (.+ b c))
      (.+ (cross a b)
	  (cross a c))))
#+nil
(c11 (vr) (vr) (vr))

(defun c12 (s a b)
  (declare (type vec a b)
	   (type num s))
  (v= (cross (.s s a) b)
      (.s s (cross a b))))

#+nil
(c12 -2.07 (v 1.7 2.6 -3.03) (v -3 2 0.2))

(defun c14 (a b c)
  (declare (type vec a b c))
  (v= (cross a (cross b c))
      (.- (.s (dot a c) b)
	  (.s (dot a b) c))))
#+nil
(c14 (vr) (vr) (vr))

(time
 (dotimes (i 300900)
   (let ((a (vr))
	 (b (vr))
	 (c (vr))
	 (s (- (random 12s0) 6s0)))
     (progn (c10 a b)
	    (c11 a b c)
	    (c12 s a b)
	    (c14 a b c)))))

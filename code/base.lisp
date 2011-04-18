(in-package :base)

(declaim (inline v copy-vec .+ .- .* ./ dot norm normalize .s))

(deftype num ()
  `single-float)

(defconstant +pi+ #.(coerce pi 'num))
(defconstant +2*pi+ #.(coerce (* 2 pi) 'num))

(defconstant +pi/2+ #.(coerce (* .5 pi) 'num))


(deftype vec ()
  `(simple-array num 1))

(deftype mat ()
  `(simple-array num 2))

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(eval-when (:compile-toplevel)
 (defun num (x)
   (declare (type number x))
   (coerce x 'single-float)))

(defmacro vec (&rest rest)
  (let ((a (gensym)))
   `(let ((,a (make-array ,(list-length rest)
			 :element-type 'num)))
      ,@(let ((i 0))
	     (declare (type (integer 0 #.(1- (expt 2 16))) i))
	     (loop for e in rest collect
		  (prog1
		       `(setf (aref ,a ,i)
			     ,(typecase e ;; FIXME once-only e?
					(fixnum (* 1s0 e))
					(num e)
					(t `(* 1s0 ,e))))
		     (incf i))))
      ,a)))

#+nil 
(vec 2 2 1)
#+nil
(vec .2 .3 1)
(defun v (&optional
	  (x 0s0)
	  (y 0s0)
	  (z 0s0))
  (the vec (vec x y z)))


(declaim (ftype (function (vec) (values vec &optional)) copy-vec))
(defun copy-vec (a)
  (let* ((n (length a))
	 (b (make-array n
			:element-type (array-element-type a))))
    (dotimes (i n)
      (setf (aref b i) (aref a i)))
    b))
#+nil
(copy-vec (v))

(progn
  (declaim (ftype (function (vec &rest t) (values vec &optional))
		  .+ .- .* ./))
  (defun .+ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (type vec e))
	(dotimes (i (length r))
	  (incf (aref r i) (aref e i))))
      r))
  (defun .- (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (type vec e))
	(dotimes (i (length r))
	  (decf (aref r i) (aref e i))))
     r))
  (defun .* (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (type vec e))
	(dotimes (i (length r))
	 (setf (aref r i) (* (aref r i) (aref e i)))))
     r))
  (defun ./ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (type vec e))
	(dotimes (i (length r))
	  (setf (aref r i) (/ (aref r i) (aref e i)))))
      r)))
#+nil
(.+ (vec 1 2 3) (vec 3 32 2) (vec 32 4 2))

(defun dot (a b)
  (declare (type vec a b))
  (let ((r 0s0))
    (declare (type num r))
    (dotimes (i (length a))
      (incf r (* (aref a i) (aref b i))))
    (the num r)))
#+nil
(dot (v 1 3 0) (v 3))

(defun norm (v)
  (declare (type vec v))
  (let ((l2 (dot v v)))
    (declare (type (single-float 0s0) l2)) ;; FIXME: write num here
    (the num (sqrt l2))))
#+nil
(norm (v 1 1 0))

(declaim (ftype (function (num vec) (values vec &optional)) .s))
(defun .s (s a)
  "scalar multiplication"
  (let ((r (v)))
    (dotimes (i (length a))
      (setf (aref r i) (* s (aref a i))))
    r))
#+nil
(.s .3 (v 1 1))


(defun normalize (v)
  (declare (type vec v))
  (let ((b (copy-vec v)))
    (the vec (.s (/ (norm v)) b))))

(defun cross (a b)
  (declare (type vec a b))
  (the vec
    (v (- (* (vy a) (vz b))
	  (* (vz a) (vy b)))
       (- (* (vz a) (vx b))
	  (* (vx a) (vz b)))
       (- (* (vx a) (vy b))
	  (* (vy a) (vx b))))))

#+nil
(cross (v 1 0 0) (v 0 1 0))


(declaim (inline vx vy vz))
(defun vx (v)
  (declare (type vec v))
  (the num (aref v 0)))

(defun vy (v)
  (declare (type vec v))
  (the num (aref v 1)))

(defun vz (v)
  (declare (type vec v))
  (the num (aref v 2)))

(defun v-spherical (theta phi)
  "Convert spherical coordinates into cartesian."
  (declare (type num theta phi) 
   #+nil (type (single-float 0s0 #.(/ (coerce pi 'num) 4)) theta)
	   #+nil (type (single-float #.(coerce (- pi) 'num) #.(coerce pi 'num)) phi))
  (let* ((st (sin theta)))
    (the vec
      (v (* st (cos phi))
	 (* st (sin phi))
	 (cos theta)))))


;; (defun check-unit-vector (&rest rest)
;;   (declare (ignore rest)))

;; (defun check-range (min max &rest rest)
;;   (declare (ignore min max rest)))


(defun check-unit-vector (&rest rest)
  ;; turn off for faster speed
  (dolist (e rest)
    (unless (< (abs (- (norm e) 1)) 1s-6)
      (error "vector isn't normalized"))))

#+nil
(check-unit-vector (v 1 0 0))

(defun check-range (min max &rest rest)
#+nil  (declare (type num min max))
  (dolist (e rest)
    (declare (type num e))
    (unless (< min e max)
      (error "range check failed"))))

(defun req (&optional name)
  (error "Required argument ~@[~S~] missing" name))

(defun m (a b c d e f g h i)
  (declare (type num a b c d e f g h i))
  (the mat
    (make-array '(3 3)
		:element-type 'num
		:initial-contents (list (list a b c) (list d e f) (list g h i)))))

(defun rotation-matrix (angle vect)
  "Create matrix that rotates by ANGLE radians around the direction
 VECT. VECT must be normalized."
  (declare (type num #+nil(single-float 0s0 #.+2*pi+) angle)
           (type vec vect))
  (check-unit-vector vect)
  (let* ((u (vx vect))
         (v (vy vect))
         (w (vz vect))
         (c (cos angle))
         (s (sin angle))
         (1-c (- 1 c))
         (su (* s u))
         (sv (* s v))
         (sw (* s w)))
    (the mat
      (m (+ c (* 1-c u u))
	 (+ (* 1-c u v) sw)
	 (- (* 1-c u w) sv)
	 
	 (- (* 1-c u v) sw)
	 (+ c (* 1-c v v))
	 (+ (* 1-c v w) su)
	 
	 (+ (* 1-c u w) sv)
	 (- (* 1-c v w) su)
	 (+ c (* 1-c w w))))))

(defun chop (a)
  (let* ((r (make-array (array-dimensions a) 
			:element-type (array-element-type a)))
	 (a1 (sb-ext:array-storage-vector a))
	 (r1 (sb-ext:array-storage-vector r)))
    (dotimes (i (length r1))
      (when 
	  (setf (aref r1 i) (if (< (abs (aref a1 i)) 1s-7)
				0s0
				(aref a1 i)))))
    r))

#+nil
(chop
 (rotation-matrix (/ +pi+ 2) (v 0 0 1)))

(defun determinant (m)
  (declare (type mat m))
  (with-arrays (m)
   (let ((a (m 0 0))
	 (b (m 0 1))
	 (c (m 0 2))
	 (d (m 1 0))
	 (e (m 1 1))
	 (f (m 1 2))
	 (g (m 2 0))
	 (h (m 2 1))
	 (i (m 2 2)))
     (the num
       (+ (* a (- (* e i) (* f h)))
	  (* b (- (* f g) (* d i)))
	  (* c (- (* d h) (* e g))))))))

(defun transpose (a)
  (declare (type mat a))
  (let ((res (m 0 0 0 0 0 0 0 0 0)))
    (loop for i below 3 do
          (loop for j below 3 do
                (setf (aref res i j)
                      (aref a j i))))
    (the mat res)))

(defun m* (matrix vect)
  "Multiply MATRIX with VECT."
  (declare (type mat matrix)
           (type vec vect))
  (let ((res (v)))
    (dotimes (i 3)
          (dotimes (j 3)
	    (incf (aref res i) 
		  (* (aref matrix i j) (aref vect j)))))
    (the vec res)))

;; rotate e_x around z, result should be -e_y=(0 -1 0)

#+nil
(chop
 (m*
  (rotation-matrix (/ +pi+ 2) (v 0 0 1))
  (v 1 0 0)))


#+nil
(let ((n 13))
 (loop for i below n collect
      (m* (rotation-matrix (* i (/ +2*pi+ n)) (normalize (v 1))) (v 0 0 1))))


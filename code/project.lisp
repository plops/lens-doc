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

(with-asy "/dev/shm/projection-test.asy"
  ;; draw rays starting from a billboard in the out-of-focus nucleus
  ;; through a target point and through coverslip
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  
  (let ((rays 90)
	(ne 1.3s0))
    (multiple-value-bind (x1 y1 hx hy s c r)
	(prepare-out-of-focus (v 0 0 .1) (v 2 0 2) :ne ne) 
      ;; green nucleus
      (asy "draw(shift(~a)*scale3(~f)*unitsphere,green);" (coord s) r)
      ;; lines on surface of cone
      (dotimes (i rays)
	(let* ((e (calc-e (* (/ +2*PI+ rays) i) x1 y1 hx hy s))
	       (f (.- c e))
	       (hf (normalize f))
	       (normal (v 0 0 1))
	       (i (intersect-plane e hf (v) normal))
	       (ht (refract-plane hf normal (/ ne 1.5))))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord e)
	       (coord i)
	       (coord (.+ i (.s 5s0 ht)))))))))

;; normal is directed towards sample (opposite i)
(defun refract-objective-detection (p i n f center normal)
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
(refract-objective-detection (v .1 0 -3.1) (v 0 0 1) 1.5s0 2s0 (v) (v 0 0 -1))


(with-asy "/dev/shm/objective.asy"
  ;; draw rays starting from point p through the objective into the
  ;; back focal plane (plot with asymptote -V /dev/shm/objective.asy)
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let ((rays 13)
	(p (v .1 0 -3.1)))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (r e) (refract-objective-detection p dir
								1.5 2.0 (v) (v 0 0 -1))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord p)
	       (coord e)
	       (coord (.+ e (.s 1s0 r)))))))))

(defun refract-objective-illumination (start dir n f center normal)
  ;; normal is directed in the direction of i
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

(with-asy "/dev/shm/objective-compare.asy"
  ;; check that raytracing from the back focal plane into the sample
  ;; gives the same rays as the reverse
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,-3.2)--(0,0,1),blue);")
  (let* ((rays 31)
	 (n 1.8)
	 (f 2.0)
	 (start (v .3 0 (+ -.2 (- (* n f)))))
	 (obj-c (v))
	 (obj-n (v 0 0 -1))
	 (bfp-c (v 0 0 f))
	 (bfp-n obj-n))
    (asy "draw(shift((0,0,~a))*rotate(90,(1,0,0))*scale3(~a)*unitcircle3);"
	 (* -1 n f) (* n f))
    (line (v 0 0 (- f)) (v 2 0 (- f)))
   (dotimes (i rays)
      (let* ((dir (v-spherical (* 2.0s0 (/ (+ .000001 (- i (floor rays 2))) rays))
			       0s0)))
	(multiple-value-bind (r e)
	    (refract-objective-detection start dir n f obj-c obj-n)
	  (let ((dy (v 0 .01 0))
		(bfp (intersect-plane e r bfp-c bfp-n)))
	    (line-colored "red" start e 
			  (.- e (.s 4s0 dy))
			  (.+ e r) bfp)
	    (line-colored "cyan" start (.+ start (v 0 0 (- (* n f) f)))
			  (v))
	    (multiple-value-bind (rr ee)
		(refract-objective-illumination bfp (.s -1s0 r)
						n f obj-c obj-n)
	      (line-colored "blue" 
			    bfp ee (.+ ee rr)))
	   ))))))



(defun refract-thin-lens (start i f center normal)
  ;; normal is directed in the opposite direction of i
  (declare (type vec start i normal center)
	   (type num f))
  (let* ((lens-hit (intersect-plane start i center normal))
	 (rho (.- lens-hit center))
	 (cos-theta (dot i normal))
	 (r (.- (.s (/ f cos-theta) i)
		rho)))
    (values (normalize r) lens-hit)))

#+nil
(refract-thin-lens (v -.1 0 -2) (normalize (v .1 0 1))
		   2.0 (v) (v 0 0 1))

(with-asy "/dev/shm/tubelens.asy"
  ;; draw rays starting from point p and going through a thin lens
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let ((rays 13)
	(p (v -.1 0 -2)))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (r e) (refract-thin-lens p dir 2.0 (v) (v 0 0 1))
	  (asy "draw((~a)--(~a)--(~a));" 
	       (coord p)
	       (coord e)
	       (coord (.+ e (.s 1s0 r)))))))))


(with-asy "/dev/shm/microscope.asy"
  ;; draw rays along sample->objective->tubelens->camera using
  ;; asymptote
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let* ((rays 4)
	 (p (v .1 0 -3.1))
	 (obj-c (v))
	 (obj-n (v 0 0 -1))
	 (obj-f 2.0)
	 (tube-f 10.0)
	 (tube-c (v 0 0 (+ obj-f tube-f)))
	 (tube-n (v 0 0 1)))
    (dotimes (i rays)
      (dotimes (j rays)
       (let* ((dir (normalize (v (/ (- i (floor rays 2)) rays)
				 (/ (- j (floor rays 2)) rays)
				 1))))
	 (multiple-value-bind (r e) (refract-objective-detection p dir
								 1.5 obj-f obj-c obj-n)
	   (multiple-value-bind (rr ee) (refract-thin-lens e r tube-f tube-c tube-n)
	     (asy "draw((~a)--(~a)--(~a)--(~a));" 
		  (coord p)
		  (coord e)
		  (coord ee)
		  (coord (.+ ee (.s (* 1.2 tube-f) rr)))
		  ))))))))

;; normal opposing dir
(defun aberrate-index-plane (start dir center normal ne/n)
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

(with-asy "/dev/shm/objective-aberrate.asy"
  ;; draw rays from sample (water) through coverslip and objective
  ;; into back focal plane
  (asy "import three;")
  (asy "size(1000,1000);")
  ;; coordinate axes
  (asy "draw((0,0,0)--(1,0,0),red);")
  (asy "draw((0,0,0)--(0,1,0),green);")
  (asy "draw((0,0,0)--(0,0,1),blue);")
  (let* ((rays 13)
	 (p (v .1 0 -3.1))
	 (obj-center (v))
	 (obj-normal (v 0 0 -1))
	 (obj-f 2.0)
	 (n 1.5)
	 (ne 1.3)
	 (slide-center (v 0 0 (* -.3 n obj-f)))
	 (slide-normal obj-normal))
    (dotimes (i rays)
      (let* ((dir (normalize (v (- (/ i rays) .5) 0 1))))
	(multiple-value-bind (dir! start! f!)
	    (aberrate-index-plane p dir slide-center slide-normal (/ ne n))
	  (multiple-value-bind (r e) 
	      (refract-objective-detection start! dir!
					   n obj-f obj-center obj-normal)
	   #+nil(line start!
		  (.+ start! (.s 2s0 dir!)))
	    (line p f! e (.+ e r))))))))

(loop for embedding-depth in '(.001 .003 .1 .03 .01) do
     ;; draw rays through a microscope
     ;; sample->coverslip->objective->tubelens->camera as asy file and graph
     ;; focus displacement in dependence of back focal plane position
     ;; with gnuplot
     (with-asy "/dev/shm/microscope-aberrate.asy"
       (asy "import three;")
       (asy "size(3000,3000);")
       (asy "currentprojection=orthographic((0,10,0),(1,0,0));")
       ;; coordinate axes
       (asy "draw((0,0,0)--(1,0,0),red);")
       (asy "draw((0,0,0)--(0,1,0),green);")
       (asy "draw((0,0,0)--(0,0,1),blue);")
       (let* ((rays 30)
	      (obj-c (v))
	      (obj-n (v 0 0 -1))
	      (obj-f 2.61)
	      (tube-f 16.0)
	      (tube-c (v 0 0 (+ obj-f tube-f)))
	      (tube-n (v 0 0 1))
	      (n 1.52)
	      (ne 1.33)
					; (embedding-depth .001) 
	      (displaced-focus (v 0 0 (+ (* -1 ne embedding-depth)
					 (- (* n (- obj-f embedding-depth))))))
	      (slide-center (.+ displaced-focus
				(.s (* -1 n embedding-depth) obj-n)))
	      (slide-normal obj-n)
	      
	      (start (.+ displaced-focus (v 0 0 0))))
	 (line-colored "blue" slide-center ;; interface between embedding and immersion 
		       (.+ slide-center (v 2 0 0)))
	 (line (v 0 0 (* -1 n obj-f)) ;; focus for ne=n
	       (v 1 0 (* -1 n obj-f)))
	 (line-colored "red" displaced-focus ;; this is the focus for an on-axis ray 
		       (.+ displaced-focus (v -2 0 0)))
	 (asy "draw(shift((0,0,~a))*rotate(90,(1,0,0))*scale3(~a)*unitcircle3);"
	      (* -1 n obj-f) (* n obj-f)) ;; gaussian sphere
	 (line (.- tube-c (v 3 0 0))      ;; tubelens
	       (.+ tube-c (v 3 0 0)))
	 
	 (let ((cam (v 0 0 (+ obj-f tube-f tube-f)))) ;; image
	   (line (.- cam (v 1 0 0)) 
		 (.+ cam (v 1 0 0)))
	   (line-colored "red" (v 0 0 -4) cam))
	 (with-open-file (gp-zoomed "/dev/shm/focus-displacement-zoom.gp"
				    :direction :output
				    :if-exists :supersede
				    :if-does-not-exist :create)
	   (with-open-file (gp "/dev/shm/focus-displacement.gp"
			       :direction :output
			       :if-exists :supersede
			       :if-does-not-exist :create)
	     (format gp "set terminal pdf; set output \"/dev/shm/focus-displacement.pdf\"; set grid;
set xlabel \"bfp ray intersection/mm\";
set ylabel \"focus displacement/mm\";
set title \"Water depths 1, 3, 10, 30, 100 um\";
plot \"/dev/shm/focus-displacement.dat\" u 1:2 w l;")
	     (format gp-zoomed "set terminal pdf; set output \"/dev/shm/focus-displacement-zoomed.pdf\";
 set grid;
set xlabel \"bfp ray intersection/mm\";
set ylabel \"focus displacement/mm\";
set yrange [-1:1];
set title \"Water depths 1, 3, 10, 30, 100 um\";
plot \"/dev/shm/focus-displacement.dat\" u 1:2 w l;"))) 
	 ;; rm focus-displacement.dat; for i  in  focus-displacement_*.dat; do cat $i >> focus-displacement.dat; echo >> focus-displacement.dat ; done; for i in *.gp ;do gnuplot $i ;done

	 (with-open-file (splot (format nil "/dev/shm/focus-displacement_~a.dat" embedding-depth)
				:direction :output
				:if-exists :supersede :if-does-not-exist :create)
	   (dotimes (i rays)
	     (let* ((dir (v-spherical (* 1.5 (atan (* 3.1s0 (/ (- i (floor rays 2)) rays))))
				      0s0)))
	       (multiple-value-bind (dir! start! f!)
		   (aberrate-index-plane start dir slide-center slide-normal (/ ne n))
		 (multiple-value-bind (r e)
		     (refract-objective-detection start! dir!
						  n obj-f obj-c obj-n)
		   (multiple-value-bind (rr ee) 
		       (refract-thin-lens e r tube-f tube-c tube-n)
		     (line start f! e ee (.+ ee (.s (* 1.2 tube-f) rr)))
		     (let ((bfp (intersect-plane e r (v 0 0 obj-f) (v 0 0 1)))
			   (hit (intersect-plane ee rr (v 0 0 (+ obj-f tube-f tube-f)) (v 0 0 1))))
		       (format splot "~a ~a~%" (vx bfp) (vx hit))))))))))))
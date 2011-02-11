;; evaluate bead sectioning
;; look at the stack
;; and compare images that were illuminated from different angles
(eval-when (:compile-toplevel :load-toplevel :execute)
  (push #P"../../0102/woropt-cyb-0628/" asdf:*central-registry*)
  (require :vol)
  (push #P"../../0102/woropt-cyb-0628/" asdf:*central-registry*)
  (require :run-ics))
;(declaim (optimize (speed 3) (safety 1) (debug 1)))
(defpackage :bead-eval
  (:use :cl :vol)
  (:export
   #:extract
   #:list->stack
   #:max-intensity-projection
   #:s*
   #:flatten
   #:with-copy-and-1d
   #:clamp-a
   #:byte-swap
   #:make-gauss
   #:make-gauss3))
(in-package :bead-eval)

(defparameter *tmp* "/home/martin/tmp/a0208/")
(defparameter *data* "/home/martin/cyberpower/d0126/")
(defparameter *secdir* "6f-60x16.3ms-delay20s-section/")
(defparameter *ang* "6e-60x16.3ms-delay20s-n8/")
(defparameter *voldir* "6-dz1um/")

(defmacro dir (&rest rest)
  `(concatenate 'string *data* ,@rest))
(defmacro tmp (&rest rest)
  `(concatenate 'string *tmp* ,@rest))

(defun loadim (fn)
  (extract (vol:convert-2-ub16/df-mul (byte-swap (vol:read-pgm fn)))
	   :start '(324 486)
	   :size '(300 300)))

(defun writeim (fn im)
  (vol:write-pgm fn (vol:normalize-2-df/ub8 im)))

;;;; VOLUME OF SECTIONS
;; load all darkimages for the section volume and integrate
#+nil
(time
 (progn
   (defparameter *avg-bg*
     (let* ((fns (directory (dir *voldir* "*dark.pgm")))
	   (bgs (loop for e in fns collect (loadim e)))
	    (n (length bgs))
	    (avg-bg (make-array (array-dimensions (first bgs))
			       :element-type 'double-float)))
       (destructuring-bind (y x) (array-dimensions (first bgs))
	 (dolist (e bgs)
	  (vol:do-region ((j i) (y x))
	    (incf (aref avg-bg j i) (* (/ 1d0 n) (aref e j i)))))
	 avg-bg)))
   (vol:write-pgm (tmp "0bg-avg.pgm") (vol:normalize-2-df/ub8 *avg-bg*))))



(defparameter *gauss* nil)
;; prepare gaussian for smoothing
#+nil
(progn
  (setf *gauss* (make-gauss *avg-bg*))
  (vol:write-pgm (tmp "0gauss.pgm")
		 (vol:normalize-2-cdf/ub8-realpart *gauss*)))

(defun blur (img &optional (kernel *gauss*))
  (vol:convert-2-csf/sf-realpart
   (vol:convolve-circ-2-csf 
    kernel
    (vol:convert-2-sf/csf-mul img))))

;; smooth background image
#+nil
(vol:write-pgm (tmp "0bg-filt.pgm") 
	       (vol:normalize-2-df/ub8
		(blur *avg-bg*)))

;; load sections (front of the array doesn't contain anything)
#+nil
(defparameter *e*
  (let* ((fns (subseq (directory (dir *voldir* "*sec*.pgm")) 13))	 
	 (stackl (loop for e in fns collect (loadim e)))
	 (stack (list->stack stackl)))
    stack))


;; construct gaussian 3GAU
#+nil
(progn
  (defparameter *g3* (make-gauss3 *e* :sigma-x-pixel  5d0))
  (save-stack-ub8 (tmp "g3") (normalize-3-cdf/ub8-realpart *g3*)))

;; convolve sectioned stack
#+nil
(progn
 (defparameter *ge*
   (convert-3-cdf/df-realpart
    (convolve-circ-3-cdf *g3* (convert-3-df/cdf-mul *e*))))
 (save-stack-ub8 (tmp "ge")
		(normalize-3-df/ub8 *ge*)))

#+nil
(progn
 (write-pgm (tmp "0xz.pgm")
	    (normalize-2-df/ub8
	     (resample-2-df 
	      (cross-section-xz *e* 212)
	      1s0 1s0 1s0 (/ 1s0 20))))
 (write-pgm (tmp "0xz-filt.pgm")
	    (normalize-2-df/ub8
	     (resample-2-df 
	      (cross-section-xz *ge* 212)
	      1s0 1s0 1s0 (/ 1s0 20))))
 (write-pgm (tmp "0xz-g.pgm")
	    (normalize-2-df/ub8
	     (resample-2-df 
	      (cross-section-xz (convert-3-cdf/df-realpart *g3*))
	      1s0 1s0 1s0 (/ 1s0 20)))))



;; do maximum intensity projection of raw data
#+nil
(progn
  (defparameter *me* (max-intensity-projection *e*))
  (vol:write-pgm (tmp "0max-proj.pgm") 
		 (vol:normalize-2-df/ub8 
		  (vol::.- *me* *avg-bg*))))

;; maximum intensity projection of filtered data
#+nil
(progn
  (defparameter *meg* (max-intensity-projection *ge*))
  (vol:write-pgm (tmp "0max-proj-filt.pgm") (vol:normalize-2-df/ub8 
					     *meg*)))

;(load "/home/martin/floh/0102/woropt-cyb-0628/run-ics.lisp")

#+nil
(format t "~a~%"(run-ics::biggest-part 
  (run-ics::point-list-sort (run-ics::nuclear-seeds *ge*))
  .746))
#+nil
(save-stack-ub8 (tmp "seeds")
		(normalize-3-df/ub8
		 (run-ics::mark-nuclear-seeds *ge* :threshold .746)))

;; background for long integration of section
#+nil
(progn
  (defparameter *longint-bg*
    (loadim (dir *secdir* "000-0-dark.pgm")))
  (vol:write-pgm (tmp "1sec-bg.pgm")
		 (vol:normalize-2-df/ub8 *longint-bg*)))

;; store section minus darkimage
#+nil
(progn
  (defparameter *longint-section* 
    (.- (loadim (dir *secdir* "000-2-section.pgm"))
	*longint-bg*))
  (write-pgm (tmp "2sec-bg.pgm")
	     (normalize-2-df/ub8 *longint-section*)))

;; store widefield minus darkimage
#+nil
(progn
  (defparameter *longint-widefield*
    (vol::.- (loadim (dir *secdir* "000-0-bright.pgm"))
	     *longint-bg*))
  (write-pgm (tmp "2wide.pgm")
		 (vol:normalize-2-df/ub8 *longint-widefield*)))

;; blur section and dark image and check the height at the in-focus bead center 
#+nil
(let* ((a (blur *longint-widefield*))
       (b (blur *longint-section*))
       (x 170)
       (y 211)
       (eta (/ (aref a y x)
	     (aref b y x)))
       (scale 13d0)
       (c (.- a (s* scale b))))
					;(setf (aref c y x) .0d0)
  (format t "~a~%" eta)
  ;; the sectioned image is nearly 12 times darker
  ;; as expected. I can't see why eta should be negative, though.
  (write-pgm (tmp "3blur-wf.pgm") (normalize-2-df/ub8 a))
  (write-pgm (tmp "3blur-sec.pgm") (normalize-2-df/ub8 b))
  (write-pgm (tmp "2subtract-blur.pgm")
	     (normalize-2-df/ub8 c))
  (write-pgm (tmp "2subtract.pgm")
	     (normalize-2-df/ub8 (.- *longint-widefield* 
				     (s* scale *longint-section*)))))


#+nil
(progn
  (defparameter *ang-bg* (loadim (dir *ang* "snap065.pgm")))
  (writeim (tmp "3bg.pgm") *ang-bg*))
#+nil
(progn
  (defparameter *ang-bright* (loadim (dir *ang* "snap064.pgm")))
  (writeim (tmp "3bright.pgm") *ang-bright*))

(defparameter *out-of-focus* nil)

;; scale all images so that the bead has the same intensity
;; integrate over region, leaving out the in focus bead
#+nil
(time
 (destructuring-bind (h w) (array-dimensions *longint-section*)
  (let* ((gauss (make-gauss *longint-section* 4d0)) 
	 (bs (s* 1.4d9 (blur *longint-section* gauss)))
	 (n 8)
	 (mosaic (make-array (list (* n h) (* n w)) :element-type 'double-float)))
   (setf *out-of-focus* nil)
    (dotimes (i 64)
      (let* ((fn (format nil "snap~3,'0d.pgm" i)) 
	     (img (.- (loadim (dir *ang* fn))
		      *ang-bg*))
	     (y 211)
	     (x 170)
	     (radius 17d0)
	     (radius2 (* radius radius))
	     (bimg (blur img gauss))
	     (eta (/ (aref *ang-bright* y x)
		   (aref bimg y x)))
	     #+nil (sub (.- bimg bs)))
	(let ((f (flatten bimg)))
	#+nil  (format t "~a~%" (list 'i i
				 'min-b (reduce #'min f)
				 'max-b (reduce #'max f)
				 'eta eta
				 'bs/bimg (/ (aref bs y x)
					     (aref bimg y x)))))
	(when (<= 2 i)
	 (let ((oj (floor i n))
	       (oi (mod i n))
	       (sum 0d0))
	   (do-region ((v u) (h w))
	     (let* ((xx (- u x))
		    (yy (- v y))
		    (r2 (+ (* 1d0 xx xx) (* 1d0 yy yy)))
		    (val (if (< radius2 r2)
			     (* eta (aref bimg v u))
			     0d0)))
	       (setf (aref mosaic
			   (+ (* h oj) v)
			 (+ (* w oi) u))
		   val)
	       (incf sum val)))
	   (format t "~a~%" sum)
	   (push (list sum (list oj oi))
		 *out-of-focus*)))))
    (vol:write-pgm (tmp "3mosaic.pgm")
		       (normalize-2-df/ub8 mosaic)))))

(defun scale-point-list (ls)
 (let* ((so (sort (copy-seq ls) #'< :key #'first))
       (mi (first (first so)))
       (ma (first (first (last so)))))
   (list ma mi)))

(scale-point-list *out-of-focus*)

#+nil
(reduce #'max (flatten *section*))

#+nil ;; this creates quite a lot useless and wrong code
(defmacro do-reg (vol &body body)
  `(etypecase ,vol
     ((simple-array * 1) (destructuring-bind (x) (array-dimensions ,vol)
			   (dotimes (i x)
			     ,@body)))
     ((simple-array * 2) (destructuring-bind (y x) (array-dimensions ,vol)
			   (dotimes (j y)
			     (dotimes (i x)
			       ,@body))))
     ((simple-array * 3) (destructuring-bind (z y x) (array-dimensions ,vol)
			   (dotimes (k z)
			     (dotimes (j y)
			       (dotimes (i x)
				 ,@body)))))))
(defun list->stack (ls)
  (declare (values (simple-array double-float 3) &optional))
  (format t "~a~%" 'list->stack)
  (let* ((n (length ls))
	 (f (first ls))
	 (dims (array-dimensions f))
	 (stack (make-array (push n dims)
			    :element-type (array-element-type f)))
	 (k 0))
    (destructuring-bind (y x) (array-dimensions f)
      (dolist (e ls)
	(format t "~a~%" k)
	(vol:do-region ((j i) (y x))
	  (setf (aref stack k j i) (aref e j i)))
	(incf k)))
    stack))


(defun extract (a &key (start nil) (center nil) (size nil))
  (let ((b (make-array size :element-type (array-element-type a))))
    (destructuring-bind (yy xx) size
      (let* ((centert (if center
			 center
			 (mapcar #'(lambda (x) (floor x 2))
				 (array-dimensions a)))) 
	    (startt (if start
			start
			(mapcar #'- centert
				(mapcar #'(lambda (x) (floor x 2))
					size)))))
	(destructuring-bind (sy sx) startt
	  (vol:do-region ((j i) (yy xx))
	    (setf (aref b j i) (aref a (+ sy j) (+ sx i)) ))
	  b)))))

(defun max-intensity-projection (vol)
  (declare ((simple-array double-float 3) vol)
	   (values (simple-array double-float 2) &optional))
  (destructuring-bind (z y x) (array-dimensions vol)
   (let* ((m (make-array (list y x) :element-type 'double-float)))
     (vol:do-region ((j i) (y x))
       (let ((ma (aref vol (1- z) j i)))
	 (dotimes (k (1- z))
	   (let ((v (aref vol k j i)))
	     (when (< ma v)
	       (setf ma v))))
	 (setf (aref m j i) ma)))
     m)))

(defmacro with-copy-and-1d (arrays &body body)
  "Provides an array NEW-A and displaced 1D array NEW-A1."
  ;; for each of supplied arrays, e.g. a, create new-a and new-a1. 
  (flet ((make-let (ls fmt cmd)
	   (loop for e in ls collect 
		(list (intern (format nil fmt (symbol-name e)) *package*)
		      ;; maybe use alexandria:format-symbol
		      (funcall cmd e)))))
    (let ((news (make-let
		 arrays "NEW-~a"
		 #'(lambda (x) `(make-array (array-dimensions ,x)
				       :element-type (array-element-type ,x)))))
	  (news1 (make-let
		  arrays "NEW-~a1"
		  #'(lambda (x) `(make-array (reduce #'* (array-dimensions ,x))
					:element-type (array-element-type ,x)
					:displaced-to ,(intern (format nil "NEW-~a" x)
							       *package*)))))
	  (olds1 (make-let
		  arrays "~a1"
		  #'(lambda (x) `(make-array (reduce #'* (array-dimensions ,x))
					:element-type (array-element-type ,x)
					:displaced-to ,x)))))
     `(let* (,@news
	     ,@news1
	     ,@olds1)
	,@body))))

(defun s* (scale vol &optional (offset 0d0))
  (declare ((simple-array double-float 2) vol)
	   (double-float scale offset)
	   (values (simple-array double-float 2) &optional))
  (with-copy-and-1d (vol)
    (declare (type (simple-array double-float 2) vol new-vol)
	     (type (array double-float 1) vol1 new-vol1))
    (dotimes (i (length vol1))
      (setf (aref new-vol1 i) (* scale (- (aref vol1 i)
					  offset))))
    new-vol))

(defun flatten (a)
  (declare (array a)
	   (values array &optional))
  (with-copy-and-1d (a)
    (dotimes (i (length a1))
      (setf (aref new-a1 i) (aref a1 i)))
    new-a1))

(defun clamp-a (a &key (scale 1d0) (offset 0d0))
  (declare ((simple-array double-float 2) a)
	   (double-float scale offset)
	   (values (array double-float 2) &optional))
  (with-copy-and-1d (a)
    (declare ((simple-array double-float 2) a new-a)
	     ((array double-float 1) a1 new-a1))
    (dotimes (i (length a1))
      (setf (aref new-a1 i) (let ((v (* scale (- (aref a1 i) offset))))
			      (cond ((< v 0d0) 0d0)
				    ((< 255d0 v) 255d0)
				    (t v)))))
    new-a))

;; for some reason the camera delivers the wrong endianess
(defun byte-swap (a)
  (declare ((simple-array (unsigned-byte 16) 2) a)
	   (values (simple-array (unsigned-byte 16) 2) &optional))
  (with-copy-and-1d (a)
    (declare ((array (unsigned-byte 16) 1) a1 new-a1))
      (dotimes (i (length a1))
	(let ((v (aref a1 i)))
	  (setf (aref new-a1 i) (+ (* 255
				      (ldb (byte 8 0) v))
				   (ldb (byte 8 8) v)))))
      new-a))

(defun make-gauss (a &optional (sigma-pixel 3d0))
 (let* ((m (make-array (array-dimensions a)
		       :element-type '(complex double-float)))
	(sigma (/ sigma-pixel (array-dimension a 1)))
	(arg (/ (* -2 sigma sigma))))
   (destructuring-bind (y x) (array-dimensions m)
    (do-region ((j i) (y x))
      (let* ((ii (* (/ 1d0 x) (- i (floor x 2))))
	     (jj (* (/ 1d0 y) (- j (floor y 2))))
	     (r2 (+ (* ii ii) (* jj jj))))
	(setf (aref m j i) (complex (exp (* arg r2)))))))
   m))

(declaim (inline sq))
(defun sq (x)
  (declare (single-float x)
	   (values single-float &optional))
  (* x x))

(defun make-gauss3 (a &key 
		    (sigma-x-pixel 3s0)
		    (sigma-y-pixel sigma-x-pixel)
		    (sigma-z-pixel (/ sigma-x-pixel 20s0)))
  (declare (type (simple-array single-float 3) a)
	   (type single-float sigma-x-pixel sigma-y-pixel sigma-z-pixel))
  (destructuring-bind (z y x) (array-dimensions a)
   (let* ((m (make-array (array-dimensions a)
			 :element-type '(complex single-float)))
	  (sx (/ (* (sqrt 2s0) sigma-x-pixel)))
	  (sy (/ (* (sqrt 2s0) sigma-y-pixel)))
	  (sz (/ (* (sqrt 2s0) sigma-z-pixel))))
     (do-region ((k j i) (z y x))
       (let* ((ii (* sx (- i (floor x 2))))
	      (jj (* sy (- j (floor y 2))))
	      (kk (* sz (- k (floor z 2))))
	      (r2 (+ (sq ii) (sq jj) (sq kk))))
	 (setf (aref m k j i) (complex (- (exp (- r2))
					  (exp (* -1.4s0 r2)))))))
     (the (simple-array (complex single-float) 3) m))))


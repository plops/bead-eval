;; evaluate bead sectioning
;; look at the stack
;; and compare images that were illuminated from different angles
(push #P"/home/martin/floh/0102/woropt-cyb-0628/" asdf:*central-registry*)
#.(require :vol)
(declaim (optimize (speed 3) (safety 1) (debug 1)))
(defpackage :bead-eval
  (:use :cl))
(in-package :bead-eval)

(defparameter *tmp* "/home/martin/tmp/a0127/")
(defparameter *data* "/home/martin/cyberpower/d0126/")
(defparameter *secdir* "6f-60x16.3ms-delay20s-section/")
(defparameter *longint* "6e-60x16.3ms-delay20s-n8/")
(defparameter *voldir* "6-dz1um/")

(defmacro dir (&rest rest)
  `(concatenate 'string *data* ,@rest))
(defmacro tmp (&rest rest)
  `(concatenate 'string *tmp* ,@rest))

(defun loadim (fn)
  (vol:convert-2-ub16/df-mul (byte-swap (vol:read-pgm fn))))

(defun writeim (fn im)
  (vol:write-pgm fn (vol:normalize-2-df/ub8 im)))


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
   (vol:write-pgm (tmp "bg-avg.pgm") (vol:normalize-2-df/ub8 *avg-bg*))))


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
#+nil
(progn
  (sb-ext:gc :full t)
  (defparameter *gauss* (let* ((m (make-array (array-dimensions *avg-bg*)
					     :element-type '(complex double-float)))
			      (sigma .01d0)
			      (arg (/ (* -2 sigma sigma))))
			 (do-reg m
			   (let* ((ii (* (/ 1d0 x) (- i (floor x 2))))
				  (jj (* (/ 1d0 y) (- j (floor y 2))))
				  (r2 (+ (* ii ii) (* jj jj))))
			     (setf (aref m j i) (complex (exp (* arg r2))))))
			 m))
  (vol:write-pgm (tmp "gauss.pgm")
		 (vol:normalize-2-cdf/ub8-realpart *gauss*)))


#+nil
(vol:write-pgm (tmp "bg-filt.pgm") 
	       (vol:normalize-2-cdf/ub8-realpart
		(vol:convolve-circ-2-cdf 
		 *gauss*
		 (vol:convert-2-df/cdf-mul *avg-bg*))))

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

#+nil
(defparameter *e*
  (let* ((fns (directory (dir *voldir* "*sec*.pgm")))
	 (stackl (loop for e in fns collect (loadim e)))
	 (stack (list->stack stackl)))
    stack))

(defparameter *st-fn*
 (directory (dir *voldir* "*sec*.pgm")))

(defparameter *m*
 (loadim (elt *st-fn* 16)))

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

(defparameter *fb*
 (extract *m* :start '(324 486)
	  :size '(300 300)))
(vol:write-pgm (tmp "fb.pgm") (vol:normalize-2-df/ub8 *fb*))

(defparameter *fb*
 (vol:extract-bbox *m* (vol:make-bbox :start (vol::v 486d0 324.0)
				      :end (vol::v (+ 486d0 300d0)
						  (+ 324d0 300d0)))))

(defparameter *me* (max-intensity-projection *e*))
(vol:write-pgm (tmp "max-proj.pgm") (vol:normalize-2-df/ub8 (vol::.- *me* *avg-bg*)))

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

#+nil
(time
 (progn
   (sb-ext:gc :full t)
  (defparameter *maxproj*
    (let* ((fns (subseq (directory (dir *voldir* "*sec*.pgm")) 28))
	   (stackl (loop for e in fns collect (loadim e)))
	   (stack (list->stack stackl))
	   (m (make-array (array-dimensions (first stackl))
			  :element-type 'double-float)))
      (destructuring-bind (z y x) (array-dimensions stack)
	(vol:do-region ((j i) (y x))
	  (format t "~a~%" (list j i))
	  (setf (aref m j i) 
		(loop for k below z maximize (aref stack k j i))))
	m)))
  (vol:write-pgm (tmp "maxproj.pgm") 
		 (vol:normalize-2-df/ub8 (vol:.- *maxproj* *avg-bg*)))))


#+nil
(defparameter *section*
  (loadim (dir *secdir* "000-2-section.pgm")))

#+nil
(vol:write-pgm "/dev/shm/o.pgm" 
	       (vol:normalize-2-df/ub8 *section*))
#+nil
(writeim "/dev/shm/bg.pgm" (loadim (dir *longint* "snap062.pgm")))
#+nil
(time
 (let ((bg (loadim (dir *longint* "snap065.pgm"))))
   (dotimes (i 64)
     (let* ((fn (format nil "snap~3,'0d.pgm" i)) 
	    (img (loadim (dir *longint* fn))))
       (vol:write-pgm (concatenate 'string "/dev/shm/2" fn)
		      (clamp (vol:.- img (s* .7d0 *section*))
			     :scale 1d0 :offset 1000d0))
       (reduce #'max (flatten img))))))
#+nil
(reduce #'max (flatten *section*))

(defun s* (scale vol &optional (offset 0d0))
  (declare ((simple-array double-float 2) vol)
	   (double-float scale offset)
	   (values (simple-array double-float 2) &optional))
  (with-copy-and-1d (vol)
    (declare ((simple-array double-float 2) vol new-vol)
	     ((simple-array double-float 1) vol1 new-vol1))
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

(defun clamp (a &key (scale 1d0) (offset 0d0))
  (declare ((simple-array double-float 2) a)
	   (double-float scale offset)
	   (values (array double-float 2) &optional))
  (with-copy-and-1d (a)
    (declare ((simple-array double-float 2) a new-a)
	     ((simple-array double-float 1) a1 new-a1))
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

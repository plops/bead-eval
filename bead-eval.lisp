;; evaluate bead sectioning
;; look at the stack
;; and compare images that were illuminated from different angles
(push #P"/home/martin/floh/0102/woropt-cyb-0628/" asdf:*central-registry*)
(require :vol)

(defpackage :bead-eval
  (:use :cl))
(in-package :bead-eval)

(defparameter *data* "/home/martin/cyberpower/d0126/")
(defparameter *secdir* "6f-60x16.3ms-delay20s-section/")
(defparameter *longint* "6e-60x16.3ms-delay20s-n8/")
(defparameter *longint* "6-dz1um/")

(defmacro dir (&rest rest)
  `(concatenate 'string *data* ,@rest))

(defun loadim (fn)
  (vol:convert-2-ub16/df-mul (byte-swap (vol:read-pgm fn))))

(defun writeim (fn im)
  (vol:write-pgm fn (vol:normalize-2-df/ub8 im)))
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

(defun s* (vol scale &optional (offset 0d0))
  (declare ((simple-array double-float *) vol)
	   (double-float scale offset)
	   (values (simple-array double-float *) &optional))
  (with-copy-and-1d (vol)
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
	   (values (simple-array double-float 2) &optional))
  (with-copy-and-1d (a)
    (dotimes (i (length a1))
      (setf (aref new-a1 i) (let ((v (* scale (- (aref a1 i) offset))))
			      (cond ((< v 0d0) 0d0)
				    ((< 255d0 v) 255d0)
				    (t v)))))
    new-a))

;; for some reason the camera delivers the wrong endianess
(defun byte-swap (a)
  (declare ((simple-array (unsigned-byte 16) *) a))
  (with-copy-and-1d (a)
      (dotimes (i (length a1))
	(let ((v (aref a1 i)))
	  (setf (aref new-a1 i) (+ (* 255 (ldb (byte 8 0) v))
				   (ldb (byte 8 8) v)))))
      new-a))


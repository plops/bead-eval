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
(defmacro dir (&rest rest)
  `(concatenate 'string *data* ,@rest))
(defun loadim (fn)
  (vol:convert-2-ub16/df-mul (byte-swap (vol:read-pgm fn))))
(defun writeim (fn im)
  (vol:write-pgm fn (vol:normalize-2-df/ub8 im)))
(defparameter *section*
  (loadim (dir *secdir* "000-2-section.pgm")))


(vol:write-pgm "/dev/shm/o.pgm" 
	       (vol:normalize-2-df/ub8 *section*))
(writeim "/dev/shm/bg.pgm" (loadim (dir *longint* "snap062.pgm")))

(defparameter *blub* (loadim (dir *longint* "snap032.pgm")))

(time
 (let ((bg (loadim (dir *longint* "snap065.pgm"))))
   (dotimes (i 64)
     (let* ((fn (format nil "snap~3,'0d.pgm" i)) 
	    (img (loadim (dir *longint* fn))))
       (vol:write-pgm (concatenate 'string "/dev/shm/2" fn)
		      (clamp (vol:.- img (s* .7d0 *section*))
			     :scale 1d0 :offset 1000d0))
       (reduce #'max (flatten img))))))

(reduce #'max (flatten *section*))

(defun s* (s a)
  (declare ((simple-array double-float *) vol)
	   (double-float s)
	   (values (simple-array double-float *) &optional))
  (let* ((b (make-array (array-dimensions a)
			:element-type (array-element-type a)))
	 (b1 (make-array (reduce #'* (array-dimensions b))
			 :element-type (array-element-type b)
			 :displaced-to b))
	 (a1 (make-array (reduce #'* (array-dimensions a))
			 :element-type (array-element-type a)
			 :displaced-to a)))
    (dotimes (i (length b1))
      (setf (aref b1 i) (* s (aref a1 i))))
    b))

(defun flatten (a)
  (declare (array a)
	   (values array &optional))
  (let* ((b1 (make-array (reduce #'* (array-dimensions a))
			 :element-type (array-element-type a)))
	 (a1 (make-array (reduce #'* (array-dimensions a))
			 :element-type (array-element-type a)
			 :displaced-to a)))
   (dotimes (i (length a1))
     (setf (aref b1 i) (aref a1 i)))
   b1))

(defun clamp (a &key (scale 1d0) (offset 0d0))
  (declare ((simple-array double-float 2) a)
	   (values (simple-array (unsigned-byte 8) 2) &optional))
 #+nil (with-copy-and-1d (a) ;; idea for useful macro
    (dotimes (i (length a1))
      (setf (aref new-a1 i) (aref a1 i)))
    new-a1)
 (let* ((b (make-array (array-dimensions a)
		       :element-type '(unsigned-byte 8)))
	(b1 (make-array (reduce #'* (array-dimensions b))
			:element-type (array-element-type b)
			:displaced-to b))
	(a1 (make-array (reduce #'* (array-dimensions a))
			:element-type (array-element-type a)
			:displaced-to a)))
   (dotimes (i (length a1))
     (setf (aref b1 i) (let ((v (* scale (- (aref a1 i) offset))))
			 (cond ((< v 0d0) 0)
			       ((< 255d0 v) 255)
			       (t (floor v))))))
   b))

;; for some reason the camera delivers the wrong endianess
(defun byte-swap (a)
  (declare ((simple-array (unsigned-byte 16) *) a))
  (let* ((b (make-array (array-dimensions a)
		       :element-type (array-element-type a)))
	(b1 (make-array (reduce #'* (array-dimensions b))
			:element-type (array-element-type b)
			:displaced-to b))
	(a1 (make-array (reduce #'* (array-dimensions a))
			:element-type (array-element-type a)
			:displaced-to a)))
    (dotimes (i (length a1))
      (let ((v (aref a1 i)))
       (setf (aref b1 i) (+ (* 255 (ldb (byte 8 0) v))
			    (ldb (byte 8 8) v)))))
    b))


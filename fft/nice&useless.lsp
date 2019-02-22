"
Маленько кладбище красивого, но бесполезного кода
"

(defun v (fun fst &rest oth)
  "
  :sign: v (fun fst &rest oth)
  :desc: Applies fun elementwise to vectors.
    Silently assumed to have the same type.
    Otherwise, type is set by the first arg.
    Signals an error if given zero arguments.
  :ex  : (v #'+ #(1 2) #(3 4))
  "
  (apply #'map (type-of fst) fun fst oth))

(defmacro a.v (a fst &rest oth)
  "
  :sign: a.v (a fst &rest oth)
  :desc: multiplies vectors by given number
  :ex  : (a.v 7.0 #(1 2) #(3 4))
  "
  `(v
     #'(lambda (&rest x) (mapcar #'(lambda (x) (* ,a x)) x))
     ,fst ,@oth))

(defun nmapk (f x)
  "
  :sign: nmapk (f x)
  :args:
    :: f : f (index object) => object
    :: x : object
  :desc: Non-consing map with known index.
  :exmp: (nmapk #'(lambda (i x) (* 2 i x)) #(1 2))
  "
  ; тут не получается с let, теряется возможность применения setf
  (dotimes (i (length x))
    (setf (elt x i) (funcall f i (elt x i)) ))
  x)

; Пока не будем считать матрицу [W], а там посмотрим

(defun quad-transform-closure (k x size &key (direct T))
  (/
    (reduce #'+ (nmapk ; тут замыкание и это меедлееено
                #'(lambda (n el) (* (W (* k n) size :direct direct) el))
                x))
    (sqrt size)))

; copy-seq не копирует всё, что хотелось бы
; кажется, я начинаю понимать, почему мало кто использует лисп для расчётов..
(defun copy-array (arr &key
                   (element-type (array-element-type arr))
                   (fill-ptr     (and (array-has-fill-pointer-p arr)
                                      (fill-pointer arr)))
                   (adjustable   (adjustable-array-p arr)))
  (let* ((dims    (array-dimensions arr))
         (new-arr (make-array dims
                              :element-type element-type
                              :adjustable adjustable
                              :fill-pointer fill-ptr)))
    (dotimes (i (array-total-size arr))
      (setf (row-major-aref new-arr i) (row-major-aref arr i)))
    new-arr))

(declaim (ftype (function ((simple-vector ctype)) copy-cmlx-vector)))
(defun copy-cmlx-vector (arr size)
  (let* ((new-arr (make-array size
                              :element-type 'ctype
                              :adjustable adjustable
                              :fill-pointer fill-ptr)))
    (dotimes (i size)
      (declare (type i uint_32))
      (setf (aref new-arr i) (aref arr i)))
    new-arr))

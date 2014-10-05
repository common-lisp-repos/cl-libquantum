;;; Libquantum porting
;;; F. Saporito
;;; 25 sept 2014

(defconstant *epsilon* 0.000001)
(defconstant *pi* (* 4.0 (atan 1)))
(defconstant *num-regs* 4)
(defconstant +MAX-UNSIGNED+ 64)
(defconstant +OBJCODE_PAGE+ 65536)
(defconstant +OBJBUF_SIZE 80)
;;; objcode elements
(defconstant +INIT+ #x00)
(defconstant +CNOT+ #x01)
(defconstant +TOFFOLI+ #x02)
(defconstant +SIGMA_X+ #x03)
(defconstant +SIGMA_Y+ #x04)
(defconstant +SIGMA_Z+ #x05)
(defconstant +HADAMARD+ #x06)
(defconstant +ROT_X+ #x07)
(defconstant +ROT_Y+ #x08)
(defconstant +ROT_Z+ #x09)
(defconstant +PHASE_KICK+ #x0A)
(defconstant +PHASE_SCALE+ #x0B)
(defconstant +COND_PHASE+ #x0C)
(defconstant +CPHASE_KICK+ #x0D)
(defconstant +SWAPLEADS+ #x0E)
;;;
;;;
(defconstant +MEASURE+ #x80)
(defconstant +BMEASURE+ #x81)
(defconstant +BMEASURE_P+ #x82)
(defconstant +NOP+ #xFF)
;;;
(defglobal opstatus 0)
(defglobal objcode (make-array 255 
                               :element-type 'integer
                               :adjustable 't
                               :initial-element 0))
(defglobal position 0)
(defglobal allocated 0)
(defglobal globalfile (make-array 255
                                  :element-type 'integer
                                  :adjustable 't
                                  :initial-element #xFF))
;;; end of objcode 

(defun quantum-n2char (mu buf)
"Convert a big integer to a byte array"
(let ((size 0))
  (setq size +MAX_UNSIGNED+)
  (loop for i from 0 below size do
        (setf (aref buf i) (/ mu (ash 1 (* 8 (1- (- size i))))))
        (setq mu (rem mu (ash 1 (* 8 (1- (- size i)))))))))

(defun quantum-objcode-start ()
  (setq opstatus 1)
  (setq allocated 1)
  (setq objcode (make-array +OBJCODE_PAGE+
                            :element-type integer
                            :initial-element 0)))

(defun quantum-objcode-stop ()
  (setq opstatus 0)
  (setq allocated 0)
  (setq objcode (make-array 1
                            :element-type integer
                            :initial-element 0)))

(defun quantum-objcode-put (operation &rest vals)
)

(defun bit-vector->integer (bit-vector)
"Create a positive integer from a bit-vector."
(reduce #'(lambda (first-bit second-bit)
            (+ (* first-bit 2) second-bit))
        bit-vector))

(defun integer->bit-vector (integer)
"Create a bit-vector from a positive integer."
(labels ((integer->bit-list (int &optional accum)
           (cond ((> int 0)
                  (multiple-value-bind (i r) (truncate int 2)
                    (integer->bit-list i (push r accum))))
                 ((null accum) (push 0 accum))
                 (t accum))))
  (coerce (integer->bit-list integer) 'bit-vector)))
    
(defun ul-mult (a b)
"Perform an unsigned-long multiplication (32-bit)."
(let ((b64 #*)
      (l-b64 0)
      (j 31)
      (cc 0)
      (b32 (integer->bit-vector #xFFFFFFFF))
      (k64 0))
  (setq k64 (* a b))
  (setq b64 (integer->bit-vector k64))
  (setq l-b64 (1- (length b64)))
  (if (> l-b64 31)
      (progn
        (setq cc (- l-b64 31))
        (loop for i downfrom l-b64 to cc do
              (setf (bit b32 j) (bit b64 i))
              (decf j)))
    (setq b32 b64))
  (bit-vector->integer b32)
))

(defstruct quantum-reg
  (width 0 :type integer)
  (size 0 :type integer)
  (hashw 0 :type integer)
  (amplitude (make-array * 
                         :element-type 'complex
                         :adjustable 't
                         :initial-element #c(0.0 0.0)))
  (state (make-array * 
                     :element-type 'integer
                     :adjustable 't
                     :initial-element 0))
  (hash (make-array *
                    :element-type 'integer
                    :adjustable 't
                    :initial-element 0)) )

(defstruct quantum-matrix
  (cols 0)
  (rows 0)
  (t1 (make-array * :element-type 'complex)))

(defmacro M (m x y)
  `(aref (quantum-matrix-t1 ,m) (+ ,x (* ,y (quantum-matrix-cols ,m)))))

(defun quantum-frac-approx (a b width)
"Fractional approximation of a decimal value"
(let (  (f (coerce (/ a b) 'float))
        (g (coerce (/ a b) 'float))
        (t1 0.0)
        (num2 0)
        (i 0)
        (den2 1)
        (num1 1)
        (den1 0)
        (num 1)
        (den 1) )
  (loop do
        (progn
          (setq i (floor (+ g 0.000005)))
          (setq t1 (- i 0.000005))
          (setq g (- g t1))
          (setq g (/ 1.0 g))
          (if (> (+ (* i den1) den2) (ash 1 width)) (return))
          (setq num (+ (* i num1) num2))
          (setq den (+ (* i den1) den2))
          (setq num2 num1)
          (setq den2 den1)
          (setq num1 num)
          (setq den1 den)
          ) 
        while (> (abs (- (coerce (/ num den) 'float) f)) (/ 1.0 (* 2 (ash 1 width))))
        )
  (setq a num)
  (setq b den)
  '()
  )
)

(defun quantum-getwidth (n)
"Calculates the number of qubits required to store N"
(loop for i from 1 sum (+ i 1)
      while (< (ash 1 i) n)
      finally (return i))
)

(defun quantum-inverse-mod (n c)
"Calculate the inverse modulus of N and C"
(loop for i from 1 sum (+ i 1)
      while (/= (rem (* i c) n) 1)
      finally (return i))
)

(defmacro quantum-gcd (u v)
  `(gcd ,u ,v) )

(defmacro quantum-prob-inline (z)
  `(expt (abs ,z) 2))

(defun quantum-ipow (a b)
  (progn
    (coerce a 'integer)
    (coerce b 'integer)
    (expt a b)
    )
)

(defmacro quantum-prob (a)
  `(quantum-prob-inline ,a) )

(defun quantum-cexp (phi)
  (exp (complex 0.0 phi)) )

(defmacro quantum-new-matrix (cols rows)
  `(make-quantum-matrix 
    :rows ,rows
    :cols ,cols
    :t1 (make-array (* ,rows ,cols) 
                    :element-type 'complex
                    :initial-element #c(0.0 0.0))))

(defun quantum-print-matrix (m)
  (let ( (rows (quantum-matrix-rows m))
         (cols (quantum-matrix-cols m)))
    (loop for i from 0 below rows do
          (loop for j from 0 below cols do
                (format t "~D " (M m j i)))
          (format t "~%"))))

(defun quantum-mmult (a b)
  (let* ((b-cols (quantum-matrix-cols b))
         (b-rows (quantum-matrix-rows b))
         (a-rows (quantum-matrix-rows a))
         (c (quantum-new-matrix b-cols a-rows)))
    (loop for i from 0 below b-cols do
          (loop for j from 0 below a-rows do
                (loop for k from 0 below b-rows do
                      (setf (M c i j) (+ (M c i j) (* (M a k j) (M b i k)))))))
    c))

(defun quantum-adjoint (m)
  (let ((cols (quantum-matrix-cols m))
        (tmp #c(0.0 0.0)))
    (loop for i from 0 below cols do
          (loop for j from 0 below i do
                (setq tmp (M m i j))
                (setf (M m i j) (conjugate (M m j i)))
                (setf (M m j i) (conjugate tmp))))
    ))

(defun quantum-matrix-to-qureg (m wid)
"Convert a vector to a quantum register"
(let* ((rows (quantum-matrix-rows m))
       (l-s 0)
       (l-size (loop for i from 0 below rows sum (1+ l-s)
                     while (/= (aref (quantum-matrix-t1 m) i) #c(0.0 0.0))))
       (reg (make-quantum-reg
             :hash (make-array (ash 1 (+ wid 2))
                               :adjustable 't
                               :element-type 'integer)
             :state (make-array l-size
                                :adjustable 't
                                :element-type 'integer)
             :amplitude (make-array l-size 
                                    :element-type 'complex
                                    :adjustable 't
                                    :initial-element #c(0.0 0.0))))
       (l-j 0)
       )

  (setf (quantum-reg-width reg) wid)
  (setf (quantum-reg-size reg) l-size)
  (setf (quantum-reg-hashw reg) (+ wid 2))
  (loop for i from 0 below rows do
        (if (/= (aref (quantum-matrix-t1 m) i) #c(0.0 0.0))
            (progn
              (setf (aref (quantum-reg-state reg) l-j) i)
              (setf (aref (quantum-reg-amplitude reg) l-j) (aref (quantum-matrix-t1 m) i))
              (setf l-j (1+ l-j)) )))
  reg))

(defun quantum-new-qureg (initval wid)
"Create a new quantum register from scratch.
Allocate memory for 1 base state"
(make-quantum-reg
 :width wid
 :size 1
 :hashw (+ wid 2)
 :hash (make-array (ash 1 (+ wid 2))
                   :adjustable 't
                   :element-type 'integer
                   :initial-element 0)
 :state (make-array 1 
                    :element-type 'integer
                    :adjustable 't
                    :initial-element initval)
 :amplitude (make-array 1 
                        :element-type 'complex
                        :adjustable 't
                        :initial-element #c(1.0 0.0))))

(defun quantum-new-qureg-size (n wid)
  (make-quantum-reg
   :width wid
   :size n
   :hashw 0
   :hash (make-array 1
                     :adjustable 't
                     :element-type 'integer)
   :state (make-array 1 
                      :element-type 'integer 
                      :initial-element 0)
   :amplitude (make-array n 
                          :element-type 'complex 
                          :initial-element #c(0.0 0.0))))

(defun quantum-new-qureg-sparse (n wid)
  (make-quantum-reg
   :width wid
   :size n
   :hashw 0
   :hash (make-array 1
                     :adjustable 't
                     :element-type 'integer)
   :state (make-array n 
                      :element-type 'integer
                      :adjustable 't
                      :initial-element 0)
   :amplitude (make-array n 
                          :element-type 'complex
                          :adjustable 't
                          :initial-element #c(0.0 0.0))))

(defun quantum-qureg-to-matrix (reg)
"Convert a quantum register to a column vector"
(let ((m (quantum-new-matrix 1 (ash 1 (quantum-reg-width reg)))))
  (loop for i from 0 below (quantum-reg-size reg) do
        (setf (aref (quantum-matrix-t1 m) (aref (quantum-reg-state reg) i)) (aref (quantum-reg-amplitude reg) i)))
  m))

(defun quantum-copy-qureg (src dst)
  (setq dst (make-quantum-reg
             :width (quantum-reg-width src)
             :size (quantum-reg-size src)
             :hashw (quantum-reg-hashw src)
             :state (quantum-reg-state src)
             :amplitude (quantum-reg-amplitude src)
             :hash (quantum-reg-hash src)))
)

(defun quantum-print-qureg (reg)
"Print the contents of a quantum register to stdout"
(loop for i from 0 below (quantum-reg-size reg) do
      (format t "~D|~D> (~D) (|" (aref (quantum-reg-amplitude reg) i)
              (aref (quantum-reg-state reg) i)
              (quantum-prob-inline (aref (quantum-reg-amplitude reg) i)))
      (loop for j from (1- (quantum-reg-width reg)) downto 0 do
            (if (= (mod j 4) 3) (format t " "))
            (format t "~D" (cond 
                            ((> (logand (ash 1 j) (aref (quantum-reg-state reg) i)) 0) 1)
                            (t 0)) ))
      (format t ">)~%"))
(format t "~%")
)

(defun quantum-print-expn (reg)
"Print the output of the modular exponentation algorithm"
(loop for i from 0 below (quantum-reg-size reg) do
      (format t "~D: ~D~%" i (- (aref (quantum-reg-state reg) i)
                                (* i (ash 1 (/ (quantum-reg-width reg) 2)))))
))

(defun quantum-print-hash (reg)
"Print the hash table of qureg to check if it's corrupted"
(loop for i from 0 below (ash 1 (quantum-reg-hashw reg)) do
      (if (/= i 0)
          (format t "~D: ~D ~D~%" i (1- (aref (quantum-reg-hash reg) i))
                  (aref (quantum-reg-state reg) (if (< (1- (aref (quantum-reg-hash reg) i)) 0) 0 1)
)))))

(defun quantum-add-hash (a pos reg)
  (let ((mark 0)
        (i (quantum-hash64 a (quantum-reg-hashw reg))))
    (loop while (/= (aref (quantum-reg-hash reg) i) 0) do
          (incf i)
          (if (= i (ash 1 (quantum-reg-hashw reg)))
              (if (= mark 0)
                  (progn
                    (setq i 0)
                    (setq mark 1))
                (format t "ERROR : HASH FULL~%"))))
    (setf (aref (quantum-reg-hash reg) i) (1+ pos))))

(defun quantum-reconstruct-hash (reg)
(block hashing-rec
  (if (= (quantum-reg-hashw reg) 0)
      (return-from hashing-rec))
  (loop for i from 0 below (ash 1 (quantum-reg-hashw reg)) do
        (setf (aref (quantum-reg-hash reg) i) 0))
  (loop for i from 0 below (quantum-reg-size reg) do
        (quantum-add-hash (aref (quantum-reg-state reg) i) i reg))))

(defun quantum-addscratch (bits reg)
"Add additional space to a qureg. It is initialized to zero and can be used as scratch space. Note that the space gets added at the LSB"
(let ((l 0))
  (setf (quantum-reg-width reg) (+ (quantum-reg-width reg) bits))
  (loop for i from 0 below (quantum-reg-size reg) do
        (setq l (ash (aref (quantum-reg-state reg) i) bits))
        (setf (aref (quantum-reg-state reg) i) l))
))

(defun quantum-kronecker (reg1 reg2)
"Compute the Kronecker product of two quantum registers"
(let* ((reg-width (+ (quantum-reg-width reg1) (quantum-reg-width reg2)))
       (reg-size (* (quantum-reg-size reg1) (quantum-reg-size reg2)))
       (reg-hashw (+ 2 reg-width))
       (l1 0)
       (reg (make-quantum-reg
             :width reg-width
             :size reg-size
             :hashw reg-hashw
             :state (make-array reg-size :element-type 'integer :adjustable 't)
             :amplitude (make-array reg-size 
                                    :element-type 'complex
                                    :adjustable 't
                                    :initial-element #c(0.0 0.0))
             :hash (make-array (ash 1 reg-hashw)
                               :adjustable 't
                               :element-type 'integer))))
  (loop for i from 0 below (quantum-reg-size reg1) do
        (loop for j from 0 below (quantum-reg-size reg2) do
              (setq l1 (+ j (* i (quantum-reg-size reg2))))
              (setf (aref (quantum-reg-state reg) l1)
                    (logior (ash (aref (quantum-reg-state reg1) i) (quantum-reg-width reg2))
                            (aref (quantum-reg-state reg2) j)))
              (setf (aref (quantum-reg-amplitude reg) l1)
                   (* (aref (quantum-reg-amplitude reg1) i) (aref (quantum-reg-amplitude reg2) j)))))
  reg))

(defun quantum-state-collapse (pos value reg)
"Reduce the state vector after measurement or partial trace"
(let* ((l-siz 0)
       (d 0.0)
       (lpat 0)
       (rpat 0)
       (j 0)
       (out (make-quantum-reg))
       (pos2 (ash 1 pos)))
  (loop for i from 0 below (quantum-reg-size reg) do
        (if (or (and (logand (aref (quantum-reg-state reg) i) pos2)
                     value)
                (and (lognand (aref (quantum-reg-state reg) i) pos2)
                     (lognot value)))
            (progn
              (setq d (+ d (quantum-prob-inline (aref (quantum-reg-amplitude reg) i))))
              (setq l-siz (1+ l-siz)))))
  (setf out (make-quantum-reg
             :size l-siz
             :width (- (quantum-reg-width reg) 1)
             :hashw (quantum-reg-hashw reg)
             :hash (quantum-reg-hash reg)
             :state (make-array l-siz :element-type 'integer :adjustable 't)
             :amplitude (make-array l-siz 
                                    :element-type 'complex
                                    :adjustable 't
                                    :initial-element #c(0.0 0.0))))
  (setq j 0)
  (loop for i from 0 below (quantum-reg-size reg) do
        (if (or (and (logand (aref (quantum-reg-state reg) i) pos2)
                     value)
                (and (lognand (aref (quantum-reg-state reg) i) pos2)
                     (lognot value)))
            (progn
              (setq rpat (loop for k from 0 below pos sum (+ rpat (ash 1 k))))
              (setq rpat (logand rpat (aref (quantum-reg-state reg) i)))
              (setq lpat (loop for k from (1- (* 8 *MAX-UNSIGNED*)) below pos by -1 sum (+ lpat (ash 1 k))))
              (setq lpat (logand lpat (aref (quantum-reg-state reg) i)))
              (setf (aref (quantum-reg-state out) j) (logior (ash lpat -1) rpat))
              (setf (aref (quantum-reg-amplitude out) j) (* (aref (quantum-reg-amplitude reg) i)
                                                            (/ 1 (sqrt d))))
              (setq j (1+ j)))))
out))

(defun quantum-hash64 (key width)
"Our 64-bit multiplicative hash function"
(let ((k32 0))
  (setq k32 (logxor (logand key #xFFFFFFFF) (ash key -32)))
  (setq k32 (ul-mult k32 #x9e370001))
  (setq k32 (ash k32 (- width 32)))
  k32))

(defun quantum-bitmask (a width bits)
"Return the reduced bitmask of a basis state"
(let ((mask 0))
  (loop for i from 0 below width do
        (if (/= (logand a (ash 1 (aref bits i))) 0)
            (incf mask (ash 1 i))))
  mask))

(defun quantum-get-state (a reg)
"Get the position of a given base state via the hash table"
(let ((i (quantum-hash64 a (quantum-reg-hashw reg))))
;  (break)
  (block hashing-0
    (if (= (quantum-reg-hashw reg) 0)
        (return-from hashing-0 a))
    (loop while (/= (aref (quantum-reg-hash reg) i) 0) do
          (if (= (aref (quantum-reg-state reg) (1- (aref (quantum-reg-hash reg) i))) a)
              (return-from hashing-0 (1- (aref (quantum-reg-hash reg) i))))
          (setq i (1+ i))
          (if (= i (ash 1 (quantum-reg-hashw reg)))
              (setq i 0)))
    (return-from hashing-0 -1))))

(defun quantum-dot-product (reg1 reg2)
"Compute the dot product of two quantum registers"
(let ((f #c(0.0 0.0))
      (j 0))
  (if (/= (quantum-reg-hashw reg2) 0)
      (quantum-reconstruct-hash reg2))

  (if (/= (aref (quantum-reg-state reg1) 0))
      (loop for i from 0 below (quantum-reg-size reg1) do
            (setq j (quantum-get-state (aref (quantum-reg-state reg1) i) reg2))
            (if (> j -1)
                (setf f (+ f (conjugate (* (aref (quantum-reg-amplitude reg1) i)
                                           (aref (quantum-reg-amplitude reg2) i)))))))
    (loop for i from 0 below (quantum-reg-size reg1) do
          (setq j (quantum-get-state i reg2))
          (if (> j -1)
              (setf f (+ f (conjugate (* (aref (quantum-reg-amplitude reg1) i)
                                         (aref (quantum-reg-amplitude reg2) i))))))))
  f))

(defun quantum-dot-product-noconj (reg1 reg2)
"Compute the dot product of two quantum registers without complex conjugation."
(let ((j 0)
      (f #c(0.0 0.0)))
  (if (/= (quantum-reg-hashw reg2) 0)
      (quantum-reconstruct-hash reg2))
  (if (= (quantum-reg-state reg2) 0)
      (setq f (loop for i from 0 below (quantum-reg-size reg1) sum (* (aref (quantum-reg-amplitude reg1) i)
                                                                      (aref (quantum-reg-amplitude reg2) (aref (quantum-reg-state reg1) i)))))
    (loop for i from 0 below (quantum-reg-size reg1) do
          (setq j (quantum-get-state (aref (quantum-reg-state reg1) i reg2)))
          (if (> j -1)
              (setq f (* (aref (quantum-reg-amplitude reg1) i) (aref (quantum-reg-amplitude reg2) j))))))
  f))

(defun quantum-scalar-qureg (r reg)
"Scalar multiplication of a quantum register. This is a purely mathematical operation without any physical meaning, so only use it if you know what you are doing."
(loop for i from 0 below (quantum-reg-size reg) do
      (setf (aref (quantum-reg-amplitude reg) i) (* (aref (quantum-reg-amplitude reg) i) r))))

(defun quantum-normalize (reg)
"Normalize a quantum register"
(let ((r 0.0))
  (setq r (loop for i from 0 below (quantum-reg-size reg) 
                sum (+ r (quantum-prob (aref (quantum-reg-amplitude reg) i)))))
  (quantum-scalar-qureg (/ 1.0 (sqrt r)) reg)))

;;; QUANTUM GATES

(defun quantum-gate1 (target m reg)
"Apply the 2x2 matrix M to the target bit. M should be unitary."
(let ((iset 0)
      (addsize 0)
      (decsize 0)
      (k 0)
      (j 0)
      (t2 #c(0.0 0.0))
      (tnot #c(0.0 0.0))
      (limit 0.0)
      (done (make-array 1 :element-type 'integer)))
  (if (/= (quantum-reg-hashw reg) 0)
      (progn
        (quantum-reconstruct-hash reg)
;        (break)
        (loop for i from 0 below (quantum-reg-size reg) do
;              (break)
              (if (= (quantum-get-state (logxor (aref (quantum-reg-state reg) i) (ash 1 target)) reg) -1)
                  (setf addsize (1+ addsize))))
;        (break)
        (setf (quantum-reg-state reg) (adjust-array (quantum-reg-state reg) (+ (quantum-reg-size reg) addsize)))
        (setf (quantum-reg-amplitude reg) (adjust-array (quantum-reg-amplitude reg) (+ (quantum-reg-size reg) addsize)))
;        (break)
        (loop for i from 0 below addsize do
              (setf (aref (quantum-reg-state reg) (+ i (quantum-reg-size reg))) 0)
              (setf (aref (quantum-reg-amplitude reg) (+ i (quantum-reg-size reg))) #c(0.0 0.0)))))
  (setq done (make-array (+ addsize (quantum-reg-size reg))
                         :element-type 'integer
                         :initial-element 0))
  (setq k (quantum-reg-size reg))
  (setq limit (* *epsilon* (/ 1.0 (ash 1 (quantum-reg-width reg)))))
  (loop for i from 0 below (quantum-reg-size reg) do
        (if (= (aref done i) 0)
            (progn
              (setq iset (logand (aref (quantum-reg-state reg) i) (ash 1 target)))
              (setq tnot #c(0.0 0.0))
              (setq j (quantum-get-state (logxor (aref (quantum-reg-state reg) i) (ash 1 target)) reg))
              (if (>= j 0)
                  (setf tnot (aref (quantum-reg-amplitude reg) j)))
              (setq t2 (aref (quantum-reg-amplitude reg) i))
              (if (/= iset 0)
                  (setf (aref (quantum-reg-amplitude reg) i) 
                        (+ (* (aref (quantum-matrix-t1 m) 2) tnot)
                           (* (aref (quantum-matrix-t1 m) 3) t2)))
                (setf (aref (quantum-reg-amplitude reg) i)
                      (+ (* (aref (quantum-matrix-t1 m) 0) t2)
                         (* (aref (quantum-matrix-t1 m) 1) tnot))))
              (if (>= j 0)
                  (progn
                    (if (/= iset 0)
                        (setf (aref (quantum-reg-amplitude reg) j)
                              (+ (* (aref (quantum-matrix-t1 m) 0) tnot)
                                 (* (aref (quantum-matrix-t1 m) 1) t2)))
                      (setf (aref (quantum-reg-amplitude reg) j)
                            (+ (* (aref (quantum-matrix-t1 m) 2) t2)
                               (* (aref (quantum-matrix-t1 m) 3) tnot)))))
                (progn
                  (if (and (= (aref (quantum-matrix-t1 m) 1) 0) 
                           (if (/= iset 0) 't 'nil)) ; && (iset)
                      (return))
                  (if (and (= (aref (quantum-matrix-t1 m) 2) 0) 
                           (if (/= iset 0) 'nil 't)) ; && !(iset) 
                      (return))
                  (setf (aref (quantum-reg-state reg) k) (logxor (aref (quantum-reg-state reg) i) (ash 1 target)))
                  (if (/= iset 0)
                      (setf (aref (quantum-reg-amplitude reg) k) (* t2 (aref (quantum-matrix-t1 m) 1)))
                    (setf (aref (quantum-reg-amplitude reg) k) (* t2 (aref (quantum-matrix-t1 m) 2))))
                  (setq k (1+ k))))
              (if (>= j 0)
                  (setf (aref done j) 1)))))
  (setf (quantum-reg-size reg) (+ (quantum-reg-size reg) addsize))
  (if (/= (quantum-reg-hashw reg) 0)
      (progn
        (setq j 0)
        (loop for i from 0 below (quantum-reg-size reg) do
              (if (< (quantum-prob-inline (aref (quantum-reg-amplitude reg) i)) limit)
                  (progn
                    (setq j (1+ j))
                    (setq decsize (1+ decsize)))
                (if (/= j 0)
                    (progn
                      (setf (aref (quantum-reg-state reg) (- i j)) (aref (quantum-reg-state reg) i))
                      (setf (aref (quantum-reg-amplitude reg) (- i j)) (aref (quantum-reg-amplitude reg) i)))))
              (if (/= decsize 0)
                  (progn
                    (setf (quantum-reg-size reg) (- (quantum-reg-size reg) decsize))
                    (setf (quantum-reg-amplitude reg) (adjust-array (quantum-reg-amplitude reg) (quantum-reg-size reg)))
                    (setf (quantum-reg-state reg) (adjust-array (quantum-reg-state reg) (quantum-reg-size reg))))))))
  (if (> (quantum-reg-size reg) (ash 1 (1- (quantum-reg-hashw reg))))
      (format t "Warning: inefficient hash table (size ~D vs hash ~D)~%"
              (quantum-reg-size reg) (ash 1 (quantum-reg-hashw reg))))
; NO DECOHERENCE SET
; (quantum-decohere reg)
  reg))

(defun quantum-gate2 (target1 target2 m reg)
"Apply 4x4 matrix M to the bits TARGET1 TARGET2. M should be unitary. Warning: code mostly untested."
(let ((j 0)
      (l 0)
      (addsize 0)
      (decsize 0)
      (psi_sub (make-array 4 :element-type 'complex
                           :initial-element #c(0.0 0.0)))
      (base (make-array 4 :element-type 'integer
                        :initial-element 0))
      (limit 0.0)
      (bits (make-array 2 :element-type 'integer
                        :initial-element 0))
      (done (make-array * :element-type 'integer)))
  (loop for i from 0 below (ash 1 (quantum-reg-hashw reg)) do
        (setf (aref (quantum-reg-hash reg) i) 0))

  (loop for i from 0 below (quantum-reg-size reg) do
        (quantum-add-hash (aref (quantum-reg-state reg) i) i reg))

  (loop for i from 0 below (quantum-reg-size reg) do
        (if (= (quantum-get-state (logxor (aref (quantum-reg-state reg) i) (ash 1 target1)) reg) -1)
            (incf addsize))
        (if (= (quantum-get-state (logxor (aref (quantum-reg-state reg) i) (ash 1 target2)) reg) -1)
            (incf addsize)))
  (setf (quantum-reg-state reg) (adjust-array (quantum-reg-state reg) (+ (quantum-reg-size reg) addsize)))
  (setf (quantum-reg-amplitude reg) (adjust-array (quantum-reg-amplitude reg) (+ (quantum-reg-size reg) addsize)))
  (loop for i from 0 below addsize do
        (setf (aref (quantum-reg-state reg) (+ i (quantum-reg-size reg))) 0)
        (setf (aref (quantum-reg-amplitude reg) (+ i (quantum-reg-size reg))) #c(0.0 0.0)))
  (setq done (make-array (+ addsize (quantum-reg-size reg))
                         :element-type 'integer
                         :initial-element 0))
  (setq l (quantum-reg-size reg))
  (setq limit (/ (/ 1.0 (ash 1 (quantum-reg-width reg))) 1000000))
  (setf (aref bits 0) target1)
  (setf (aref bits 1) target2)
  (loop for i from 0 below (quantum-reg-size reg) do
        (if (= (aref done i) 0)
            (progn
              (setq j (quantum-bitmask (aref (quantum-reg-state reg) i) 2 bits))
              (setf (aref base j) i)
              (setf (aref base (logxor j 1)) (quantum-get-state (logxor (aref (quantum-get-state reg) i) (ash 1 target2)) reg))
              (setf (aref base (logxor j 2)) (quantum-get-state (logxor (aref (quantum-get-state reg) i) (ash 1 target1)) reg))
              (setf (aref base (logxor j 3)) (quantum-get-state (logxor (logxor (aref (quantum-get-state reg) i)
                                                                                (ash 1 target1))
                                                                        (ash 1 target2)) reg))

              (loop for j from 0 below 4 do
                    (if (= (aref base j) -1)
                        (progn
                          (setf (aref base j) l)
                          (incf l)))
                    (setf (aref psi_sub j) (aref (quantum-reg-amplitude reg) (aref base j))))

              (loop for j from 0 below 4 do
                    (setf (aref (quantum-reg-amplitude reg) (aref base j)) #c(0.0 0.0))
                    (loop for k from 0 below 4 do
                          (setf (aref (quantum-reg-amplitude reg) (aref base j)) (+ (aref (quantum-reg-amplitude reg) (aref base j))
                                                                                    (* (M m k j) (aref psi_sub k)))))
                    (setf (aref done (aref base j)) 1)))))
  (incf (quantum-reg-size reg) addsize)

  (setq j 0)
  (loop for i from 0 below (quantum-reg-size reg) do
        (if (< (quantum-prob-inline (aref (quantum-reg-amplitude reg) i)) limit)
            (progn
              (incf j)
              (incf decsize))
          (if (/= j 0)
              (progn
                (setf (aref (quantum-reg-state reg) (- i j)) (aref (quantum-reg-state reg) i))
                (setf (aref (quantum-reg-amplitude reg) (- i j)) (aref (quantum-reg-amplitude reg) i))))))
  (if (/= decsize 0)
      (progn
        (decf (quantum-reg-size reg) decsize)
        (setf (quantum-reg-amplitude reg) (adjust-array (quantum-reg-amplitude reg) (quantum-reg-size reg)))
        (setf (quantum-reg-state reg) (adjust-array (quantum-reg-state reg) (quantum-reg-size reg)))))
; (quantum-decohere reg)
))

(defun quantum-hadamard (target reg)
"Apply a hadamard gate"
(let ((m (quantum-new-matrix 2 2)))
  (setf (aref (quantum-matrix-t1 m) 0) (sqrt (/ 1.0 2.0)))
  (setf (aref (quantum-matrix-t1 m) 1) (sqrt (/ 1.0 2.0)))
  (setf (aref (quantum-matrix-t1 m) 2) (sqrt (/ 1.0 2.0)))
  (setf (aref (quantum-matrix-t1 m) 3) (- (sqrt (/ 1.0 2.0))))
  (quantum-gate1 target m reg)))

(defun quantum-walsh (width reg)
"Apply a walsh-hadamard transform"
(loop for i from 0 below width do
      (quantum-hadamard i reg)))




using System;

namespace matrix_library
{
    class Matrix
    {
        private float[,] matrix { get; set; }
        private readonly int height;
        private readonly int width;

        /// <summary>
        /// Creates a new <c>Matrix</c> object from a two-dimensional float array.
        /// </summary>
        /// <param name="array">2-dimensional float array</param>
        public Matrix(float[,] array)
        {
            matrix = array;
            height = array.GetLength(0);
            width = array.GetLength(1);
        }

        /// <summary>
        /// Creates a new null <c>Matrix</c> with specified height and width.
        /// </summary>
        /// <param name="height"></param>
        /// <param name="width"></param>
        public Matrix(int height, int width)
        {
            matrix = new float[height, width];
            this.height = height;
            this.width = width;
        }

        /// <summary>
        /// Creates a new <c>Matrix</c> object from a string,
        /// where numbers in a row are separated by space and rows are
        /// separated by semi-colon.
        /// <example>
        /// For example:
        /// <code>"3 1 0; 2 1 1; 5 6 2"</code>
        /// </example>
        /// </summary>
        /// <param name="s">Input string</param>
        public Matrix(string s)
        {
            // input is either a two-dimensional int array or a string
            // separated by semicolon, such as: " 1 2; 3 4; 5 6 "
            matrix = ParseInput(s);
            height = matrix.GetLength(0);
            width = matrix.GetLength(1);
        }

        private float[,] ParseInput(string s)
        {
            var rows = s.Trim().Split(';');
            int height = rows.Length;
            int width = rows[0].Trim().Split(' ').Length;
            float[,] matrix = new float[height, width];

            int row_index = 0;
            foreach (string row in rows)
            {
                string[] line = row.Trim().Split(' ');
                if (line.Length != width)
                    throw new FormatException("Wrong input format.");
                for (int i = 0; i < width; i++)
                {
                    matrix[row_index, i] = float.Parse(line[i]);
                }
                row_index++;
            }
            return matrix;
        }

        public float this[int row, int column]
        {
            get => matrix[row, column];
            set => matrix[row, column] = value;
        }

        private Matrix GetRow(int index)
        {
            float[,] row = new float[1, width];
            for (int i = 0; i < width; i++)
                row[0, i] = matrix[index, i];
            return new Matrix(row);
        }

        private void SetRow(Matrix row, int index)
        {
            for (int i = 0; i < width; i++)
            {
                matrix[index, i] = row[0, i];
            }
        }

        private Matrix GetColumn(int index)
        {
            float[,] col = new float[height, 1];
            for (int i = 0; i < height; i++)
                col[i, 0] = matrix[i, index];
            return new Matrix(col);
        }

        private void SetColumn(Matrix col, int index)
        {
            for (int i = 0; i < height; i++)
            {
                matrix[i, index] = col[i, 0];
            }
        }

        private void AddRows(int indexFrom, int indexTo, float multiplier = 1)
        {
            // add multiple of one row to another
            Matrix sum = (multiplier * GetRow(indexFrom)) + GetRow(indexTo);
            SetRow(sum, indexTo);
        }

        private void SwapRows(int i, int j)
        {
            var temp = GetRow(i);
            SetRow(GetRow(j), i);
            SetRow(temp, j);
        }

        private void MultiplyRowByConst(int index, float c)
        {
            var row = GetRow(index);
            row = c * row;
            SetRow(row, index);
        }

        /// <summary>
        /// Returns a new identity matrix.
        /// </summary>
        /// <param name="rank">Dimensions of the matrix</param>
        /// <returns>Identity matrix</returns>
        public static Matrix IdentityMatrix(int rank)
        {
            var matrix = new Matrix(new float[rank, rank]);
            for (int i = 0; i < rank; i++)
                matrix[i, i] = 1;
            return matrix;
        }
        
        private static Matrix Elimination(Matrix A, bool RREF = false, bool inverse = false, bool det = false)
        {
            // Gauss (-Jordan) eliminatin. Also used to calculate inverse and determinant

            Matrix m = Copy(A);
            Matrix inv = IdentityMatrix(A.height);

            int pivotColumn = 0;
            int pivotRow = 0;
            float f;
            float d = 1;

            while (pivotColumn < m.width && pivotRow < m.height)
            {
                int firstGoodPivotRow = int.MaxValue;

                for (int i = pivotRow; i < m.height; i++)
                {
                    if (m[i, pivotColumn] != 0)
                    {
                        firstGoodPivotRow = i;
                        break;
                    }
                }

                if (firstGoodPivotRow == int.MaxValue)
                    pivotColumn++;

                else
                {
                    f = 1 / m[firstGoodPivotRow, pivotColumn];
                    m.SwapRows(pivotRow, firstGoodPivotRow);

                    if (det && pivotRow != firstGoodPivotRow) d *= -1;

                    m.MultiplyRowByConst(pivotRow, f);

                    if (det) d *= 1 / f;

                    if (inverse)
                    {
                        inv.SwapRows(pivotRow, firstGoodPivotRow);
                        if (det && pivotRow != firstGoodPivotRow) d *= -1;
                        inv.MultiplyRowByConst(pivotRow, f);
                        if (det) d *= 1 / f;
                    }

                    for (int i = 0; i < m.height; i++)
                    {
                        if (!RREF && i <= pivotRow) continue;
                        if (RREF && i == pivotRow) continue;

                        f = m[i, pivotColumn] / m[pivotRow, pivotColumn];
                        m.AddRows(pivotRow, i, -f);

                        if (inverse)
                            inv.AddRows(pivotRow, i, -f);
                    }

                    pivotRow++;
                    pivotColumn++;
                }
            }

            if (inverse)
                return inv;

            if (det)
            {
                for (int i = 0; i < m.height; i++)
                    d *= m[i, i];
                return new Matrix($"{d}");
            }

            return m;
        }

        /// <summary>
        /// Performs the Gaussian elimination on a <c>Matrix</c> object.
        /// </summary>
        /// <param name="A"></param>
        /// <returns>Returns a new <c>Matrix</c> in the row echelon form.</returns>
        public static Matrix REF(Matrix A)
        {
            return Elimination(A);
        }

        /// <summary>
        /// Performs the Gauss-Jordan elimination on a <c>Matrix</c> object.
        /// </summary>
        /// <param name="A"></param>
        /// <returns>Returns a new <c>Matrix</c> in the reduced row echelon form.</returns>
        public static Matrix RREF(Matrix A)
        {
            return Elimination(A, true);
        }

        private static bool IsZero(Matrix row)
        {
            // returns true if all values in the row are zero
            for (int i = 0; i < row.width; i++)
            {
                if (row[0, i] != 0)
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Calculates the rank of a matrix.
        /// </summary>
        /// <param name="A">Input matrix</param>
        /// <returns>Integer representing the number of linearly independent rows.</returns>
        public static int Rank(Matrix A)
        {
            var m = Elimination(A);

            for (int i = 0; i < m.height; i++)
            {
                if (IsZero(m.GetRow(i)))
                    return i;
            }
            return m.height;
        }

        /// <summary>
        /// Calculates the inverse of a <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Regular matrix</param>
        /// <returns>New inverted matrix</returns>
        public static Matrix Inverse(Matrix A)
        {
            if (A.height != A.width || Rank(A) != A.height)
                throw new Exception("Can only invert regular matrix.");
            return Elimination(A, true, true);
        }

        /// <summary>
        /// Calculates the determinant of a <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <returns>Float value - determinant of a matrix.</returns>
        public static float Det(Matrix A)
        {
            if (A.height != A.width)
                throw new Exception("Can only calculate determinant of a square matrix.");

            else return Elimination(A, false, false, true)[0, 0];
        }

        /// <summary>
        /// Calculates the trace of a <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>The sum of elements on the diagonal.</returns>
        public static float Trace(Matrix A)
        {
            float tr = 1;
            for (int i = 0; i < Math.Min(A.height, A.width); i++)
                tr += A[i, i];
            return tr;
        }

        /// <summary>
        /// Makes a new copy of a <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Matrix to be copied</param>
        /// <returns>New Matrix object.</returns>
        public static Matrix Copy(Matrix A)
        {
            // creates a copy of the matrix
            float[,] copy = new float[A.height, A.width];
            for (int i = 0; i < A.height; i++)
            {
                for (int j = 0; j < A.width; j++)
                {
                    copy[i, j] = A[i, j];
                }
            }
            return new Matrix(copy);
        }

        private static float Norm(Matrix v)
        {
            // eucleidian norm
            return (float)Math.Sqrt(ScalarProduct(v, v));
        }

        /// <summary>
        /// Calculates the eigenvalues of a <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>For matrices with full column rank returns a vector
        /// of all eigenvalues. Othrewise returns the eigenvalue with
        /// largest absolute value.</returns>
        public static Matrix EigenValues(Matrix A)
        {
            try
            {
                var eigVals = QR_Algorithm(A, 50).Diag();
                return eigVals;
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine("Couldn't calculate all eigenvalues.\nCalculating the largest one...");
                var maxEigVector = PowerIteration(A, 50);
                return new Matrix($"{(A * maxEigVector)[0, 0] / maxEigVector[0, 0]}");
            }
        }

        /// <summary>
        /// Converts a vector (<c>Matrix</c> whose one dimension is
        /// equal to one) to a <c>float[]</c>.
        /// </summary>
        /// <param name="A">Vector</param>
        /// <returns>One-dimensional array of floats</returns>
        public static float[] ToArray(Matrix A)
        {
            if (A.height != 1 && A.width != 0)
                throw new Exception("Cannot convert this Matrix to one-dimensional array.");

            if (A.width == 1)
                A = A.Transpose();

            var array = new float[Math.Max(A.width, A.height)];
            for (int i = 0; i < A.height; i++)
                array[i] = A[i, 0];

            return array;
        }

        private static Matrix PowerIteration(Matrix A, int iter)
        {
            Matrix guess = new Matrix(A.height, 1);
            Matrix next_guess = new Matrix(A.height, 1);
            guess = IdentityMatrix(A.height).GetColumn(0);

            for (int i = 0; i < iter; i++)
            {
                next_guess = (1 / Norm(A * guess)) * (A * guess);
                guess = next_guess;
            }
            return guess;
        }

        private static Matrix QR_Algorithm(Matrix A, int iter)
        {
            // algorithm for computing eigenvalues
            Matrix A_next = new Matrix(A.height, A.width);
            for (int i = 0; i < iter; i++)
            {
                var decomp = QR_Decomp(A);
                var Q = decomp[0];
                var R = decomp[1];
                A_next = R * Q;
                A = A_next;
            }
            return A_next;
        }

        private static Matrix[] QR_Decomp(Matrix A)
        {
            // A = Q*R, s.t Q orthogonal and R upper triang.
            Matrix Q = GramSchmidt(A);
            Matrix R = Q.Transpose() * A;
            return new Matrix[2] { Q, R };
        }

        /// <summary>
        /// Does Gram-Schmidt process for columns of a <c>Matrix</c> object.
        /// </summary>
        /// <param name="A">Input Matrix</param>
        /// <returns>Returns new orthogonal matrix whose columns are
        /// the orthonormal basis of S(A)</returns>
        private static Matrix GramSchmidt(Matrix A)
        {
            // A must be square and its columns linearly independent for QR decomp.
            // u_k = a_k - (sum j=1 -> (k-1) (proj(a_k, u_j)))    u_k is perpendicular to u_1..u_(k-1)
            // e_k = u_k / ||u_k||                                normalise

            if (Rank(A.Transpose()) != A.width)
                throw new Exception("Columns of input matrix must be linearly independent.");

            Matrix q = new Matrix(A.height, A.width);

            // fill up q with orthogonal basis
            var sum = new Matrix(A.height, 1);
            for (int i = 0; i < A.width; i++)
            {
                sum.ZeroOut();
                for (int j = 0; j < i; j++)
                    sum += Proj(A.GetColumn(i), q.GetColumn(j));
                var u = A.GetColumn(i) - sum;
                q.SetColumn(u, i);
            }

            // normalise columns of q
            for (int i = 0; i < A.width; i++)
            {
                var u = q.GetColumn(i);
                var e = (1 / Norm(u)) * u;
                q.SetColumn(e, i);
            }

            return q;
        }

        private static Matrix Proj(Matrix a, Matrix u)
        {
            // projection of vector a onto vector u
            // (<u, a> / <u, u>) * u
            float c = ScalarProduct(u, a) /
                      ScalarProduct(u, u);
            return c * u;
        }

        /// <summary>
        /// Vector of elements on the diagonal.
        /// </summary>
        public Matrix Diag()
        {
            if (height != width)
                throw new Exception("A must be a square matrix.");
            Matrix diag = new Matrix(height, 1);
            for (int i = 0; i < height; i++)
                diag[i, 0] = this[i, i];
            return diag;
        }

        private void ZeroOut()
        {
            // sets all values in matrix to zero
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    this[i, j] = 0;
        }

        /// <summary>
        /// Matrix transposition
        /// </summary>
        public Matrix Transpose()
        {
            var AT = new float[width, height];
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    AT[j, i] = this[i, j];
                }
            }
            return new Matrix(AT);
        }

        /// <summary>
        /// Calculates the standard scalar product (dot product) of two vectors
        /// (vector is a <c>Matrix</c>, whose one dimension is equal to one).
        /// </summary>
        /// <param name="a">vector</param>
        /// <param name="b">vector</param>
        /// <returns>Float that is equal to scalar product of the input vectors.</returns>
        public static float ScalarProduct(Matrix a, Matrix b)
        {
            // standard scalar (dot) product of two vectors

            if (a.height == 1)
                a = a.Transpose();

            if (b.height == 1)
                b = b.Transpose();

            if (a.height != b.height || a.width != 1 || b.width != 1)
                throw new Exception("Cannot do scalar product.");

            float res = 0;

            for (int i = 0; i < a.height; i++)
                res += a[i, 0] * b[i, 0];
            
            return res;
        }

        /// <summary>
        /// Converts a <c>Matrix</c> object to a string. Rounds the values
        /// to three decimal places.
        /// </summary>
        /// <returns>String that represents the table of values in the matrix</returns>
        public override string ToString()
        {
            string output = "";
            double temp;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    temp = Math.Round(matrix[i, j], 3);
                    if (temp == -0)
                        temp = 0;
                    output += temp.ToString() + " ";
                }
                output += "\n";
            }
            return output;
        }

        /// <summary>
        /// Performs addition on two <c>Matrix</c> objects with same dimensions.
        /// </summary>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            if ((A.height != B.height) || (A.width != B.width))
                throw new FormatException("Wrong matrix dimensions");

            var sum = new float[A.height, A.width];
            for (int i = 0; i < A.height; i++)
            {
                for (int j = 0; j < A.width; j++)
                {
                    sum[i, j] = A[i, j] + B[i, j];
                }
            }
            return new Matrix(sum);
        }

        /// <summary>
        /// Performs subtraction of two <c>Matrix</c> objects with same dimensions.
        /// </summary>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            return A + ((-1) * B);
        }

        /// <summary>
        /// Scalar multiple of a <c>Matrix</c>.
        /// </summary>
        /// <param name="c">constant</param>
        /// <param name="A">Matrix object</param>
        public static Matrix operator *(float c, Matrix A)
        {
            var multiple = new float[A.height, A.width];
            for (int i = 0; i < A.height; i++)
            {
                for (int j = 0; j < A.width; j++)
                {
                    multiple[i, j] = c * A[i, j];
                }
            }
            return new Matrix(multiple);
        }

        /// <summary>
        /// Matrix multiplication.
        /// Width of left <c>Matrix</c> must be equal to height of right <c>Matrix</c>.
        /// </summary>
        /// <param name="A">Left Matrix</param>
        /// <param name="B">Right Matrix</param>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.width != B.height)
                throw new Exception("Wrong matrix dimensions.");

            var product = new float[A.height, B.width];

            for (int i = 0; i < A.height; i++)
            {
                for (int j = 0; j < B.width; j++)
                {
                    product[i, j] = ScalarProduct(A.GetRow(i), B.GetColumn(j));
                }
            }
            return new Matrix(product);
        }

        static void Main(string[] args)
        {
            // examples
            var A = new Matrix("4 3 2 2; 0 1 0 -2; 1 -1 3 3; 2 3 1 1");

            Console.WriteLine("A =");
            Console.WriteLine(A);

            Console.WriteLine($"determinant = {Det(A)}\n");

            Console.WriteLine("REF:");
            Console.WriteLine(REF(A));

            Console.WriteLine("RREF:");
            Console.WriteLine(RREF(A));

            Console.WriteLine("A^(-1)");
            Console.WriteLine(Inverse(A));

            Console.WriteLine("Eigenvalues:");
            Console.WriteLine(EigenValues(A));

            Console.WriteLine("3 * A =");
            Console.WriteLine(3 * A);

            Console.WriteLine("A * A =");
            Console.WriteLine(A * A);

            Console.WriteLine("QR decomposition:\nA =\n");
            var QR = QR_Decomp(A);
            Console.WriteLine(QR[0]);
            Console.WriteLine("*\n");
            Console.WriteLine(QR[1]);

            Matrix B = new Matrix("2 5 1 4; 4 1 6 3; 5 3 7 2; 1 0 2 4");

            Console.WriteLine("B =");
            Console.WriteLine(B);

            Console.WriteLine("Matrix multiplication:\nA * B =");
            Console.WriteLine(A * B);
        }
    }
}

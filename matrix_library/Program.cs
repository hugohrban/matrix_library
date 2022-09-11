using System;

namespace matrix_library
{
    class Matrix
    {
        private float[,] matrix { get; set; }
        private readonly int height;
        private readonly int width;

        public Matrix(float[,] input)
        {
            matrix = input;
            height = input.GetLength(0);
            width = input.GetLength(1);
        }

        public Matrix(int height, int width)
        {
            matrix = new float[height, width];
            this.height = height;
            this.width = width;
        }

        public Matrix(string s)
        {
            // input is either a two-dimensional int array or a string separated by semicolon
            // in format such as: " 1 2; 3 4; 5 6 "
            matrix = ParseInput(s);
            height = matrix.GetLength(0);
            width = matrix.GetLength(1);
        }

        private float this[int row, int column]
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

        private float[,] ParseInput(string s)
        {
            // parse string input TODO put this in constructor
            // TODO "I_n" creates identity matrik of rank n

            var rows = s.Trim().Split(';');
            int height = rows.Length;
            int width = rows[0].Trim().Split(' ').Length;
            float[,] matrix = new float[height, width];

            int row_index = 0;
            foreach (string row in rows) //TODO as a normal for loop
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

        public static Matrix IdentityMatrix(int rank)
        {
            var matrix = new Matrix(new float[rank, rank]);
            for (int i = 0; i < rank; i++)
                matrix[i, i] = 1;
            return matrix;
        }

        //public static int[] Eigenvalues()
        //{
        //    
        //}

        //public static int[,] Eigenvectors()
        //{
        //    
        //}

        //public static int Rank(Matrix M)
        //{
        //   
        //}

        private static Matrix Elimination(Matrix A, bool RREF = false, bool inverse = false)
        {
            // Gauss (-Jordan) elimination. Also used to calculate inverse
            // TODO maybe change of base

            Matrix m = Copy(A);
            Matrix inv = IdentityMatrix(A.height);

            int pivotColumn = 0;
            int pivotRow = 0;
            float f;

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
                    f = 1 / m[pivotRow, pivotColumn];
                    m.SwapRows(pivotRow, firstGoodPivotRow);
                    m.MultiplyRowByConst(pivotRow, f);

                    if (inverse)
                    {
                        inv.SwapRows(pivotRow, firstGoodPivotRow);
                        inv.MultiplyRowByConst(pivotRow, f);
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
            return m;
        }

        public static int Rank(Matrix m)
        {
            var A = Elimination(m);
            bool IsZero(Matrix row)
            {
                for (int i = 0; i < A.width; i++)
                {
                    if (row[0, i] != 0)
                        return false;
                }
                return true;
            }

            for (int i = 0; i < A.height; i++)
            {
                if (IsZero(A.GetRow(i)))
                    return i;
            }
            return A.height;
        }

        public static Matrix REF(Matrix A)
        {
            return Elimination(A);
        }

        public static Matrix Inverse(Matrix A)
        {
            if (A.height != A.width || Rank(A) != A.height)
                throw new Exception("cannot invert this.");
            return Elimination(A, true, true);
        }

        public static Matrix RREF(Matrix A)
        {
            return Elimination(A, true);
        }

        private static Matrix Copy(Matrix A)
        {
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

        public static float det(Matrix A)
        {
            if (A.height != A.width)
                throw new Exception("Can only calculate determinant of a square matrix.");

            if (A.height == 2)
            {
                return (A[0, 0] * A[1, 1]) - (A[0, 1] * A[1, 0]);
            }
            else return float.PositiveInfinity;
        }

        public static float norm(Matrix v)
        {
            return (float)Math.Sqrt(ScalarMultiple(v, v.Transpose()));
        }

        public static float[] QR_Algorithm(Matrix A)
        {
            // algorithm for computing eigen(values/vectors)
            var decomp = QR_Decomp(A);
            var Q = decomp[0];
            var R = decomp[1];

            //TODO

            return null;
        }

        public static Matrix[] QR_Decomp(Matrix A)
        {
            // A = Q*R, s.t Q orthogonal and R upper triang.
            Matrix Q = GramSchmidt(A);
            Matrix R = Q.Transpose() * A;
            return new Matrix[2] { Q, R };
        }

        private static Matrix GramSchmidt(Matrix A)
        {
            // does G-S process for columns of A. Returns new orthogonal matrix which columns are the orth. basis
            Matrix q = new Matrix(A.height, A.width);

            //TODO

            // u_k = a_k - (sum j=1 -> (k-1) (proj(a_k, u_j)))    u_k is perpendicular to u_1..u_(k-1)
            // e_k = u_k / ||u_k||                                normalise

            //Matrix step(int i)
            //{
            //    Matrix sum = new Matrix(A.height, 1);
            //    for (int j = 0; j < i; j++)
            //        sum += proj(A.GetColumn(i), q.GetColumn(j));
            //    return A.GetColumn(i) - sum;
            //}

            // fill up q with orthogonal basis
            var sum = new Matrix(A.height, 1);
            for (int i = 0; i < A.width; i++)
            {
                //sum = Matrix(A.height, 1);
                sum.ZeroOut();
                for (int j = 0; j < i; j++)
                        sum += proj(A.GetColumn(i), q.GetColumn(j));
                var u = A.GetColumn(i) - sum;
                q.SetColumn(u, i);
            }

            // normalise columns of q
            for (int i = 0; i < A.width; i++)
            {
                var u = q.GetColumn(i);
                var e = (1 / norm(u)) * u;
                q.SetColumn(e, i);
            }

            return q;
        }

        private static Matrix proj(Matrix a, Matrix u)
        {
            // projection of vector a onto vector u. (<u, a> / <u, u>) * u
            float c = ScalarMultiple(u.Transpose(), a) /
                      ScalarMultiple(u.Transpose(), u);
            return c * u;
        }

        private void ZeroOut()
        {
            // sets all values in matrix to zero
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    this[i, j] = 0;
        }

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

        public static float ScalarMultiple(Matrix a, Matrix b)
        {
            // standard scalar (dot) product of two vectors

            if ((a.height != 1) || (b.width != 1) || (a.width != b.height)) 
                throw new Exception("cannot do scalar product of this");

            float res = 0;

            for (int i = 0; i < a.width; i++)
                res += a[0, i] * b[i, 0];
            
            return res;
        }

        public override string ToString()
        {
            string output = "";
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    if (matrix[i, j] == -0)
                        matrix[i, j] = 0;
                    output += matrix[i, j].ToString() + " ";
                }
                output += "\n";
            }
            return output;
        }

        public static Matrix operator +(Matrix a, Matrix b)
        {
            if ((a.height != b.height) || (a.width != b.width))
                throw new FormatException("incorrect matrix dimensions");

            var sum = new float[a.height, a.width];
            for (int i = 0; i < a.height; i++)
            {
                for (int j = 0; j < a.width; j++)
                {
                    sum[i, j] = a[i, j] + b[i, j];
                }
            }
            return new Matrix(sum);
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            return A + ((-1) * B);
        }

        public static Matrix operator *(float c, Matrix a)
        {
            var multiple = new float[a.height, a.width];
            for (int i = 0; i < a.height; i++)
            {
                for (int j = 0; j < a.width; j++)
                {
                    multiple[i, j] = c * a[i, j];
                }
            }
            return new Matrix(multiple);
        }

        public static Matrix operator *(Matrix a, Matrix b)
        {
            if (a.width != b.height)
                throw new Exception("wrong matrix dimensions");

            var product = new float[a.height, b.width];

            for (int i = 0; i < a.height; i++)
            {
                for (int j = 0; j < b.width; j++)
                {
                    product[i, j] = ScalarMultiple(a.GetRow(i), b.GetColumn(j));
                }
            }
            return new Matrix(product);
        }

        static void Main(string[] args)
        {
            var A = new Matrix("1 2; 3 4");
            //var B = new Matrix ("0 3 1 9; 1 1 -1 1; 3 11 5 35");
            var B = new Matrix("0 1; 1 0");
            var I = Inverse(A);
            Console.WriteLine(A);
            A.ZeroOut();
            Console.WriteLine(A);
            var C = new Matrix(3, 3);
            Console.WriteLine(C);
        }
    }

    class Program { }
}

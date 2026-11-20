use std::fmt;
use Vec;

#[derive(Debug, PartialEq, Clone)]
pub struct RMatrix {
    pub elem: Vec<f32>,
    pub rows: usize,
    pub cols: usize,
}

impl fmt::Display for RMatrix {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.format_matrix())
    }
}

impl RMatrix {
    pub fn format_matrix(&self) -> String {
        if self.elem.len() != self.rows * self.cols {
            return format!(
                "[invalid RMatrix: dimensions {}x{}, {} elements]",
                self.rows,
                self.cols,
                self.elem.len()
            );
        }

        if self.rows == 0 || self.cols == 0 {
            return String::new();
        }

        let mut widths = vec![0; self.cols];
        let mut formatted = Vec::with_capacity(self.elem.len());

        for (index, value) in self.elem.iter().enumerate() {
            let text = format!("{}", value);
            let column = index % self.cols;
            widths[column] = widths[column].max(text.len());
            formatted.push(text);
        }

        let mut lines = Vec::with_capacity(self.rows);
        for row in 0..self.rows {
            let mut cells = Vec::with_capacity(self.cols);
            for column in 0..self.cols {
                let index = row * self.cols + column;
                cells.push(format!(
                    "{:>width$}",
                    formatted[index],
                    width = widths[column]
                ));
            }
            lines.push(cells.join(" "));
        }

        lines.join("\n")
    }

    /**
     * Add two matrices.
     * The first matrix is mutated and contains the summation.
     */
    pub fn mut_add(&mut self, v: &RMatrix) {
        if (self.rows, self.cols) != (v.rows, v.cols) {
            println!("The matrix layouts differ !");
            return;
        }
        if self.elem.len() != self.rows * self.cols || v.elem.len() != v.rows * v.cols {
            println!("The matrix element counts are inconsistent !");
            return;
        }

        for i in 0..self.rows * self.cols {
            self.elem[i] += v.elem[i];
        }
    }

    pub fn add(self: &RMatrix, v: &RMatrix) -> Result<RMatrix, String> {
        if (self.rows, self.cols) != (v.rows, v.cols) {
            return Err(String::from("The matrix layouts differ !"));
        }

        return Ok(RMatrix {
            elem: self
                .elem
                .iter()
                .zip(v.elem.iter())
                .map(|(u, v)| u + v)
                .collect(),
            rows: self.rows,
            cols: self.cols,
        });
    }

    /**
     * Sub operation A-B.
     * The first matrix is mutated and contains the summation.
     */
    pub fn mut_sub(&mut self, v: &RMatrix) {
        if (self.rows, self.cols) != (v.rows, v.cols) {
            println!("The matrix layouts differ !");
            return;
        }
        if self.elem.len() != self.rows * self.cols || v.elem.len() != v.rows * v.cols {
            println!("The matrix element counts are inconsistent !");
            return;
        }

        for i in 0..self.rows * self.cols {
            self.elem[i] -= v.elem[i];
        }
    }

    pub fn sub(self: &RMatrix, v: &RMatrix) -> Result<RMatrix, String> {
        if (self.rows, self.cols) != (v.rows, v.cols) {
            return Err(String::from("The matrix layouts differ !"));
        }

        return Ok(RMatrix {
            elem: self
                .elem
                .iter()
                .zip(v.elem.iter())
                .map(|(u, v)| u - v)
                .collect(),
            rows: self.rows,
            cols: self.cols,
        });
    }

    pub fn mult_by_scalar(self: &RMatrix, r: f32) -> RMatrix {
        return RMatrix {
            elem: self.elem.iter().map(|x| x * r).collect(),
            rows: self.rows,
            cols: self.cols,
        };
    }

    /**
     * Dot product of two matrices
     */
    pub fn dot(&self, v: &RMatrix) -> Result<RMatrix, String> {
        if self.cols != v.rows {
            return Err(String::from(
                "The col size of matrix 1 must equal the row size of matrix 2 !",
            ));
        }

        let mut w: RMatrix = RMatrix {
            elem: vec![0.; self.rows * v.cols],
            rows: self.rows,
            cols: v.cols,
        };
        for i in 0..self.rows {
            for j in 0..v.cols {
                let mut sum = 0.0;
                for k in 0..self.cols {
                    sum += self.elem[i * self.cols + k] * v.elem[k * v.cols + j]
                }
                w.elem[i * w.cols + j] = sum;
            }
        }
        return Ok(w);
    }

    /**
     * Swap two rows j and k of the matrix (mutates the matrix)
     * j, k: indices of the rows to swap
     * temp_row: temp vector with the number of columns of the matrix
     */
    pub fn swap_rows(&mut self, j: usize, k: usize, temp_row: &mut Vec<f32>) {
        temp_row.copy_from_slice(&self.elem[j * self.cols..(j + 1) * self.cols]);
        self.elem
            .copy_within(k * self.cols..(k + 1) * self.cols, j * self.cols);
        self.elem[k * self.cols..(k + 1) * self.cols].copy_from_slice(&temp_row);
    }

    pub fn transpose(&self) -> RMatrix {
        let nelems = self.rows * self.cols;
        let mut t: RMatrix = RMatrix {
            elem: vec![0.; nelems],
            rows: self.cols,
            cols: self.rows,
        };

        for n in 0..nelems {
            let i = n / self.rows;
            let j = n % self.rows;
            t.elem[n] = self.elem[j * self.cols + i]
        }

        return t;
    }

    /*
     * LU factorization with partial pivoting
     *
     * Code adapted from C (https://en.wikipedia.org/wiki/LU_decomposition)
     * INPUT : a - square matrix having dimension n
     *         tol - small tolerance number to detect failure when the matrix is near degenerate
     * OUTPUT: a LU matrix contains a copy of both matrices L-E and U as a=(L-E)+U such that PA=LU.
     *        The permutation matrix p is not stored as a RMatrix, but as an integer vector of size n+1
     *        containing column indexes where the permutation matrix has "1". The last element p[n]=S+n,
     *        where S is the number of row exchanges needed for determinant computation, det(p)=(-1)^S
     */
    pub fn decompose_lup(&mut self, tol: f32) -> Result<(RMatrix, Vec<usize>), String> {
        if self.rows != self.cols {
            return Err(String::from("The matrix must be square !"));
        }
        if tol <= 0.0 || !tol.is_finite() {
            return Err(String::from("The tolerance must be positive !"));
        }
        if self.elem.len() != self.rows * self.cols {
            return Err(String::from("The matrix element count is inconsistent !"));
        }
        let mut a: RMatrix = self.clone();
        let n = a.rows;

        // Permutation vector
        let mut p: Vec<usize> = vec![0; n + 1];
        for i in 0..n {
            p[i] = i
        }
        p[n] = n;

        let mut temp_row: Vec<f32> = vec![0.0f32; self.cols];
        for i in 0..n {
            // for each row
            let mut max_a: f32 = 0.0f32;
            let mut imax = i;

            for k in i..n {
                let abs_a = a.elem[k * n + i].abs();
                if abs_a > max_a {
                    // find the pivot
                    max_a = abs_a;
                    imax = k;
                }
            }

            if max_a <= tol {
                return Err(String::from("matrix is degenerate !"));
            }

            if imax != i {
                //pivot p and rows of A
                p.swap(i, imax);
                a.swap_rows(i, imax, &mut temp_row);
                //count pivots starting from n (for determinant)
                p[n] += 1;
            }

            // Gaussian elimination
            for j in i + 1..n {
                a.elem[j * n + i] /= a.elem[i * n + i];

                for k in i + 1..n {
                    a.elem[j * n + k] -= a.elem[j * n + i] * a.elem[i * n + k]
                }
            }
        }

        return Ok((a, p)); //decomposition done
    }

    /* Solve A*x = B
     *
     * INPUT: a,p filled by decompose_lup, such that PA = LU in LU*x=Pb
     *         b - rhs vector, given as a column RMatrix
     * OUTPUT: x - solution vector of a*x=b, returned as a column RMatrix
     */
    pub fn solve_lup(&self, p: &Vec<usize>, b: &RMatrix) -> Result<RMatrix, String> {
        let n = self.rows;
        if n != self.cols {
            return Err(String::from("The matrix a must be square !"));
        }
        if p.len() != n + 1 {
            return Err(String::from(
                "The p vector must have the number of rows of a + 1 !",
            ));
        }
        if b.rows != n {
            return Err(String::from(
                "The b column matrix must have the number of rows of a !",
            ));
        }
        if b.cols != 1 {
            return Err(String::from("The b matrix must be a column vector !"));
        }
        if b.elem.len() != n {
            return Err(String::from("The b element count is inconsistent !"));
        }

        // Solution vector
        let mut x = vec![0.; n];

        for i in 0..n {
            x[i] = b.elem[p[i]];
            for k in 0..i {
                x[i] -= self.elem[i * n + k] * x[k];
            }
        }

        for i in (0..n).rev() {
            for k in i + 1..n {
                x[i] -= self.elem[i * n + k] * x[k];
            }
            x[i] /= self.elem[i * n + i];
        }

        return Ok(RMatrix {
            elem: x,
            rows: n,
            cols: 1,
        });
    }

    /*
     * Invert the matrix
     * INPUT: a,p filled in decompose_lup
     * OUTPUT: the inverse of the initial matrix
     */
    pub fn inverse(self: &RMatrix, p: &Vec<usize>) -> Result<RMatrix, String> {
        let n = self.rows;
        if n != self.cols {
            return Err(String::from("The matrix must be square !"));
        }
        if p.len() != n + 1 {
            return Err(String::from("The P vector must have N+1 rows !"));
        }

        // Inverse matrix
        let mut i_a: Vec<f32> = vec![0.; n * n];

        for j in 0..n {
            for i in 0..n {
                i_a[i * n + j] = if p[i] == j { 1.0 } else { 0.0 }; // P matrix

                for k in 0..i {
                    i_a[i * n + j] -= self.elem[i * n + k] * i_a[k * n + j];
                }
            }

            for i in (0..n).rev() {
                for k in i + 1..n {
                    i_a[i * n + j] -= self.elem[i * n + k] * i_a[k * n + j];
                }

                i_a[i * n + j] /= self.elem[i * n + i];
            }
        }

        return Ok(RMatrix {
            elem: i_a,
            rows: n,
            cols: n,
        });
    }

    /* INPUT: p vector filled in decompose_lup
     * OUTPUT: return the determinant of the initial matrix
     */
    pub fn determinant(&self, p: &Vec<usize>) -> f32 {
        let n = self.rows;
        if n == 0 || n != self.cols || self.elem.len() != n * n || p.len() != n + 1 {
            return f32::NAN;
        }
        let mut det = self.elem[0] as f32;
        for i in 1..n {
            det *= self.elem[i * n + i];
        }

        let swaps = p[n].saturating_sub(n);
        if swaps % 2 == 0 {
            return det;
        } else {
            return -det;
        }
    }

    /**
     * Compute the norm of a vector
     */
    pub fn norm_vector(v: &Vec<f32>) -> f32 {
        return f32::sqrt(v.iter().map(|x| x * x).sum());
    }

    pub fn norm(v: &RMatrix) -> f32 {
        return RMatrix::norm_vector(&v.elem);
    }

    /**
     * Identity matrix
     * n : size of the matrix
     * scale: multiplicative factor (put 1. for the unit identity matrix)
     */
    pub fn eye(n: usize, scale: f32) -> RMatrix {
        let mut v = vec![0.; n * n];
        for i in 0..n {
            v[i * n + i] = scale;
        }
        return RMatrix {
            elem: v,
            rows: n,
            cols: n,
        };
    }

    /**
     * Inverse iteration method for finding an eigenvalue/eigenvector pair,
     * given an initial guess.
     *
     * INPUT: y, alpha  initial guess for the eigenvector and eigenvalue
     * INPUT: epsilon precision of the result
     * OUTPUT: an eigenvalue/eigenvector pair that is the closest to y
     */
    pub fn inverse_iter(
        self: &RMatrix,
        y: &RMatrix,
        alpha: f32,
        epsilon: f32,
    ) -> Result<(u32, f32, RMatrix), String> {
        let n = self.rows;
        if n != self.cols {
            return Err(String::from("The matrix must be square !"));
        }
        if epsilon <= 0.0 || !epsilon.is_finite() {
            return Err(String::from("The precision epsilon must be positive !"));
        }
        if y.rows != n || y.cols != 1 {
            return Err(String::from("The initial guess vector must have the dimensions (rows: n, cols: 1), where n is the size of the matrix !"));
        }
        if RMatrix::norm(y) == 0.0 {
            return Err(String::from("The initial guess vector must not be zero !"));
        }

        let mut yy = y.clone();
        let alphax_i = RMatrix::eye(n, alpha);
        let mut b = self.sub(&alphax_i)?; // A - alpha*I aka "shifted of origin"
        let (lu, permut) = b.decompose_lup(epsilon)?;

        let mut iter = 0;
        let mut distance = f32::MAX;
        let mut theta = 0.0_f32;
        const MAX_ITERATIONS: u32 = 10_000;

        while distance > epsilon * f32::abs(theta) && iter < MAX_ITERATIONS {
            iter += 1;

            let norm_y = RMatrix::norm_vector(&yy.elem);
            if norm_y == 0.0 || !norm_y.is_finite() {
                return Err(String::from(
                    "The inverse iteration produced an invalid vector !",
                ));
            }
            let v = RMatrix {
                elem: yy.mult_by_scalar(1.0 / norm_y).elem,
                rows: n,
                cols: 1,
            };

            // solve (A - sigma*I) y = v, where v = y/||y||
            yy = lu.solve_lup(&permut, &v)?;
            theta = v.transpose().dot(&yy)?.elem[0];
            if theta == 0.0 || !theta.is_finite() {
                return Err(String::from(
                    "The inverse iteration produced a zero or invalid Rayleigh quotient !",
                ));
            }
            distance = RMatrix::norm(&yy.sub(&v.mult_by_scalar(theta))?);
        }
        if iter == MAX_ITERATIONS {
            return Err(String::from("The inverse iteration did not converge !"));
        }

        return Ok((iter, alpha + 1.0 / theta, yy.mult_by_scalar(1.0 / theta)));
    }
}

#[cfg(test)]
mod tests {
    use super::RMatrix;

    fn assert_matrix_close(a: &RMatrix, b: &RMatrix) {
        assert_eq!(a.rows, b.rows);
        assert_eq!(a.cols, b.cols);
        assert_eq!(a.elem.len(), b.elem.len());
        for (x, y) in a.elem.iter().zip(b.elem.iter()) {
            assert!((x - y).abs() < 1e-5, "{:?} != {:?}", a.elem, b.elem);
        }
    }

    #[test]
    fn format_matrix_aligns_columns() {
        let matrix = RMatrix {
            elem: vec![1.0, -22.5, 300.0, 4.25],
            rows: 2,
            cols: 2,
        };

        assert_eq!(matrix.format_matrix(), "  1 -22.5\n300  4.25");
        assert_eq!(format!("{}", matrix), "  1 -22.5\n300  4.25");
    }

    #[test]
    fn add_and_sub_require_matching_layouts() {
        let a = RMatrix {
            elem: vec![1.0, 2.0],
            rows: 1,
            cols: 2,
        };
        let b = RMatrix {
            elem: vec![1.0, 2.0],
            rows: 2,
            cols: 1,
        };

        assert!(a.add(&b).is_err());
        assert!(a.sub(&b).is_err());
    }

    #[test]
    fn mutating_add_and_sub_leave_the_matrix_unchanged_on_layout_error() {
        let mut a = RMatrix {
            elem: vec![1.0, 2.0],
            rows: 1,
            cols: 2,
        };
        let b = RMatrix {
            elem: vec![3.0, 4.0],
            rows: 2,
            cols: 1,
        };

        a.mut_add(&b);
        assert_matrix_close(
            &a,
            &RMatrix {
                elem: vec![1.0, 2.0],
                rows: 1,
                cols: 2,
            },
        );

        a.mut_sub(&b);
        assert_matrix_close(
            &a,
            &RMatrix {
                elem: vec![1.0, 2.0],
                rows: 1,
                cols: 2,
            },
        );
    }

    #[test]
    fn inconsistent_mutating_element_counts_do_not_panic() {
        let mut a = RMatrix {
            elem: vec![1.0],
            rows: 1,
            cols: 2,
        };
        let b = RMatrix {
            elem: vec![3.0, 4.0],
            rows: 1,
            cols: 2,
        };

        a.mut_add(&b);
        assert_eq!(a.elem, vec![1.0]);
        a.mut_sub(&b);
        assert_eq!(a.elem, vec![1.0]);
    }

    #[test]
    fn dot_product_supports_rectangular_matrices() {
        let a = RMatrix {
            elem: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            rows: 2,
            cols: 3,
        };
        let b = RMatrix {
            elem: vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
            rows: 3,
            cols: 2,
        };

        assert_matrix_close(
            &a.dot(&b).unwrap(),
            &RMatrix {
                elem: vec![58.0, 64.0, 139.0, 154.0],
                rows: 2,
                cols: 2,
            },
        );
    }

    #[test]
    fn transpose_has_correct_dimensions_and_values() {
        let a = RMatrix {
            elem: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            rows: 2,
            cols: 3,
        };

        assert_matrix_close(
            &a.transpose(),
            &RMatrix {
                elem: vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0],
                rows: 3,
                cols: 2,
            },
        );
    }

    #[test]
    fn lu_rejects_non_square_and_degenerate_matrices() {
        let mut rectangular = RMatrix {
            elem: vec![1.0, 2.0, 3.0],
            rows: 1,
            cols: 3,
        };
        assert!(rectangular.decompose_lup(1e-6).is_err());

        let mut singular = RMatrix {
            elem: vec![1.0, 2.0, 2.0, 4.0],
            rows: 2,
            cols: 2,
        };
        assert!(singular.decompose_lup(1e-6).is_err());

        let mut zero = RMatrix {
            elem: vec![0.0],
            rows: 1,
            cols: 1,
        };
        assert!(zero.decompose_lup(0.0).is_err());
    }

    #[test]
    fn determinant_sign_uses_permutation_swap_count() {
        let mut a = RMatrix {
            elem: vec![0.0, 1.0, 1.0, 0.0],
            rows: 2,
            cols: 2,
        };
        let (lu, permut) = a.decompose_lup(1e-6).unwrap();

        assert_eq!(lu.determinant(&permut), -1.0);
    }

    #[test]
    fn solve_lup_rejects_a_non_column_rhs() {
        let mut a = RMatrix {
            elem: vec![1.0, 0.0, 0.0, 1.0],
            rows: 2,
            cols: 2,
        };
        let (lu, permut) = a.decompose_lup(1e-6).unwrap();
        let rhs = RMatrix {
            elem: vec![1.0, 2.0],
            rows: 1,
            cols: 2,
        };

        assert!(lu.solve_lup(&permut, &rhs).is_err());
    }

    #[test]
    fn inverse_returns_the_identity_product() {
        let mut a = RMatrix {
            elem: vec![4.0, 7.0, 2.0, 6.0],
            rows: 2,
            cols: 2,
        };
        let (lu, permut) = a.decompose_lup(1e-6).unwrap();
        let inverse = lu.inverse(&permut).unwrap();
        let product = a.dot(&inverse).unwrap();

        assert_matrix_close(
            &product,
            &RMatrix {
                elem: vec![1.0, 0.0, 0.0, 1.0],
                rows: 2,
                cols: 2,
            },
        );
    }

    #[test]
    fn inverse_iteration_converges_for_a_diagonal_matrix() {
        let a = RMatrix {
            elem: vec![2.0, 0.0, 0.0, 1.0],
            rows: 2,
            cols: 2,
        };
        let y = RMatrix {
            elem: vec![1.0, 1.0],
            rows: 2,
            cols: 1,
        };

        let (_, eigenvalue, eigenvector) = a.inverse_iter(&y, 0.0, 1e-6).unwrap();

        assert!((eigenvalue - 1.0).abs() < 1e-4);
        assert!(eigenvector.elem[0].abs() < 1e-4);
        assert!(eigenvector.elem[1].abs() > 0.99);
    }

    #[test]
    fn inverse_iteration_rejects_invalid_inputs_and_singular_shift() {
        let a = RMatrix {
            elem: vec![2.0, 0.0, 0.0, 1.0],
            rows: 2,
            cols: 2,
        };
        let y = RMatrix {
            elem: vec![0.0, 0.0],
            rows: 2,
            cols: 1,
        };

        assert!(a.inverse_iter(&y, 0.0, 1e-6).is_err());

        let y = RMatrix {
            elem: vec![1.0, 1.0],
            rows: 2,
            cols: 1,
        };
        assert!(a.inverse_iter(&y, 0.0, 0.0).is_err());
        assert!(a.inverse_iter(&y, 2.0, 1e-6).is_err());
    }
}

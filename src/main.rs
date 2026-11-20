mod rmatrix;
use rmatrix::RMatrix as Matrix;

fn unit_tests() {
    let mut s1 = Matrix {
        elem: vec![1., -2., 2., -1.],
        rows: 2,
        cols: 2,
    };
    let mut s2 = Matrix {
        elem: vec![3., 4., -1., 2.],
        rows: 2,
        cols: 2,
    };
    println!("s1 =\n{}", s1);
    println!("s2 =\n{}", s2);

    let mut s3 = s1.add(&s2);
    match s3 {
        Ok(s3) => {
            println!("s1+s2 =\n{}", s3.format_matrix());
            assert!(
                s3 == Matrix {
                    elem: vec![4., 2., 1., 1.],
                    rows: 2,
                    cols: 2
                }
            )
        }
        Err(e) => panic!("{}", e),
    }

    s1 = Matrix {
        elem: vec![1., -2., 2., -1.],
        rows: 2,
        cols: 2,
    };
    s2 = Matrix {
        elem: vec![3., 4., -1., 2.],
        rows: 2,
        cols: 2,
    };
    s3 = s1.sub(&s2);
    match s3 {
        Ok(s3) => {
            println!("s1-s2 =\n{}", s3.format_matrix());
            assert!(
                s3 == Matrix {
                    elem: vec![-2.0, -6.0, 3.0, -3.0],
                    rows: 2,
                    cols: 2
                }
            )
        }
        Err(e) => panic!("{}", e),
    }

    s1 = Matrix {
        elem: vec![1., -2., 2., -1.],
        rows: 2,
        cols: 2,
    };
    s2 = Matrix {
        elem: vec![3., 4., -1., 2.],
        rows: 2,
        cols: 2,
    };
    s1.mut_add(&s2);
    assert!(
        s1 == Matrix {
            elem: vec![4., 2., 1., 1.],
            rows: 2,
            cols: 2
        }
    );
    println!("s1+s2 =\n{}", s1.format_matrix());

    s1 = Matrix {
        elem: vec![1., -2., 2., -1.],
        rows: 2,
        cols: 2,
    };
    s2 = Matrix {
        elem: vec![3., 4., -1., 2.],
        rows: 2,
        cols: 2,
    };
    s1.mut_sub(&s2);
    println!("s1-s2 =\n{}", s1.format_matrix());
    assert!(
        s1 == Matrix {
            elem: vec![-2.0, -6.0, 3.0, -3.0],
            rows: 2,
            cols: 2
        }
    );

    let s4 = Matrix {
        elem: vec![1., 3., 5., 2., 4., 6.],
        rows: 2,
        cols: 3,
    };
    let s5 = Matrix {
        elem: vec![-1., 1., 2., 3., 0., -2.],
        rows: 3,
        cols: 2,
    };
    println!("s4 =\n{}", s4.format_matrix());
    println!("s5 =\n{}", s5.format_matrix());

    let mut s6 = s4.dot(&s5);
    match s6 {
        Ok(s6) => {
            println!("s4.s5 =\n{}", s6.format_matrix());
            assert_eq!(
                s6,
                Matrix {
                    elem: vec![5.0, 0.0, 6.0, 2.0],
                    rows: 2,
                    cols: 2
                }
            );
        }
        Err(e) => panic!("{}", e),
    }
    s6 = s5.dot(&s4);
    match s6 {
        Ok(s6) => {
            println!("s5.s4 =\n{}", s6.format_matrix());
            assert_eq!(
                s6,
                Matrix {
                    elem: vec![1.0, 1.0, 1.0, 8.0, 18.0, 28.0, -4.0, -8.0, -12.0],
                    rows: 3,
                    cols: 3
                }
            )
        }
        Err(e) => panic!("{}", e),
    }

    println!("s4' =\n{}", s4.transpose().format_matrix());
    println!("s5' =\n{}", s5.transpose().format_matrix());

    // LU decomposition, inversion, equation solving, determinant
    let mut square = Matrix {
        elem: vec![10., -9., 18., 7., -12., 11., -10., 10., 3.],
        rows: 3,
        cols: 3,
    };
    println!("A =\n{}", square.format_matrix());
    let (plu, permut) = square.decompose_lup(0.0001).unwrap();
    assert_eq!(
        plu,
        Matrix {
            elem: vec![
                10.0,
                -9.0,
                18.0,
                0.7,
                -5.7000003,
                -1.5999994,
                -1.0,
                -0.17543858,
                20.7193
            ],
            rows: 3,
            cols: 3
        }
    );
    println!("LU =\n{}", plu.format_matrix());
    println!("P = {:?}", permut);

    let inverse = plu.inverse(&permut);
    match inverse {
        Ok(inverse) => {
            println!("inv A =\n{}", inverse.format_matrix());
            assert_eq!(
                inverse,
                Matrix {
                    elem: vec![
                        0.12362405,
                        -0.17527518,
                        -0.09906857,
                        0.11092294,
                        -0.1778154,
                        -0.013547834,
                        0.042337,
                        0.0084673995,
                        0.04826418
                    ],
                    rows: 3,
                    cols: 3
                }
            )
        }
        Err(e) => panic!("{}", e),
    }

    let det = plu.determinant(&permut);
    assert_eq!(det, -1181.0001);
    println!("det A = {:?}", det);

    println!("Solving AX = [1,2,3]");
    let b = Matrix {
        elem: vec![1., 2., 3.],
        rows: 3,
        cols: 1,
    };
    let x = plu.solve_lup(&permut, &b).unwrap();
    assert_eq!(
        x,
        Matrix {
            elem: vec![-0.524132, -0.28535134, 0.20406434],
            rows: 3,
            cols: 1
        }
    );
    println!("X =\n{}", x.format_matrix());
    let verif = square.dot(&x);
    match verif {
        Ok(verif) => {
            assert_eq!(
                verif,
                Matrix {
                    elem: vec![1.0, 1.9999998, 3.0],
                    rows: 3,
                    cols: 1
                }
            );
            println!("Verification: A.x =\n{}", verif.format_matrix());
        }
        Err(e) => panic!("{}", e),
    }

    let eye_a = Matrix::eye(3, 3.14);
    assert_eq!(
        eye_a,
        Matrix {
            elem: vec![3.14, 0., 0., 0., 3.14, 0., 0., 0., 3.14],
            rows: 3,
            cols: 3
        }
    );

    println!("Computing eigenvalue/eigenvector for = [[1,2,0],[-1,0,1],[1,-2,-1]]");
    let c = Matrix {
        elem: vec![1., 2., 0., -1., 0., 1., 1., -2., -1.],
        rows: 3,
        cols: 3,
    };
    let y = Matrix {
        elem: vec![1., 0., 1.],
        rows: 3,
        cols: 1,
    };
    let alpha = 1.0;
    let (iter, eigenvalue, eigenvector) = c.inverse_iter(&y, alpha, 1e-7).unwrap();
    println!("iteration {:?}:", iter);
    println!("eigenvector: v =\n{}", eigenvector.format_matrix());
    println!("eigenvalue: lambda = {:?}", eigenvalue);
}

fn main() {
    unit_tests();
}

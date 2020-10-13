mod rmatrix;
use rmatrix::RMatrix as Matrix;

fn main() {

    let mut s1 = Matrix{elem: vec![1.,-2.,2.,-1.], rows: 2, cols: 2};
    let mut s2 = Matrix{elem: vec![3.,4.,-1.,2.], rows: 2, cols: 2};
    println!("s1 = {:?}", s1);
    println!("s2 = {:?}", s2);

    let mut s3 = s1.add(&s2);
    match s3{
        Ok(s3) => {
            println!("s1+s2 = {:?}", s3);
            assert!(s3 == Matrix{elem: vec![4.,2.,1.,1.], rows: 2, cols: 2})
        },
        Err(e) => panic!(e)
    }

    s1 = Matrix{elem: vec![1.,-2.,2.,-1.], rows: 2, cols: 2};
    s2 = Matrix{elem: vec![3.,4.,-1.,2.], rows: 2, cols: 2};
    s3 = s1.sub(&s2);
    match s3{
        Ok(s3) => {
            println!("s1-s2 = {:?}", s3);
            assert!(s3 == Matrix{elem: vec![-2.0, -6.0, 3.0, -3.0], rows: 2, cols: 2})
        },
        Err(e) => panic!(e)
    }

    s1 = Matrix{elem: vec![1.,-2.,2.,-1.], rows: 2, cols: 2};
    s2 = Matrix{elem: vec![3.,4.,-1.,2.], rows: 2, cols: 2};
    s1.mut_add(&s2);
    assert!(s1 == Matrix{elem: vec![4.,2.,1.,1.], rows: 2, cols: 2});
    println!("s1+s2 = {:?}", s1);

    s1 = Matrix{elem: vec![1.,-2.,2.,-1.], rows: 2, cols: 2};
    s2 = Matrix{elem: vec![3.,4.,-1.,2.], rows: 2, cols: 2};
    s1.mut_sub(&s2);
    println!("s1-s2 = {:?}", s1);
    assert!(s1 == Matrix{elem: vec![-2.0, -6.0, 3.0, -3.0], rows: 2, cols: 2});

    let s4 = Matrix{elem: vec![1.,3.,5.,2.,4.,6.], rows: 2, cols: 3};
    let s5 = Matrix{elem: vec![-1.,1.,2.,3.,0.,-2.], rows: 3, cols: 2};
    println!("s4 = {:?}", s4);
    println!("s5 = {:?}", s5);

    let mut s6 = s4.dot(&s5);
    match s6{
        Ok(s6) => {
            println!("s4.s5 = {:?}", s6);
            assert_eq!(s6, Matrix{ elem: vec![5.0, 0.0, 6.0, 2.0], rows: 2, cols: 2 });
        },
        Err(e) => panic!(e)
    }
    s6 = s5.dot(&s4);
    match s6{
        Ok(s6) => {
            println!("s5.s4 = {:?}", s6);
            assert_eq!(s6, 
                Matrix { elem: vec![1.0, 1.0, 1.0, 8.0, 18.0, 28.0, -4.0, -8.0, -12.0], 
                rows: 3, 
                cols: 3 })
        },
        Err(e) => panic!(e)
    }

    println!("s4' = {:?}", s4.transpose());
    println!("s5' = {:?}", s5.transpose());

    // LU decomposition, inversion, solving, determinant
    let mut square = Matrix{
        elem: vec![10.,-9.,18.,7.,-12., 11.,-10.,10.,3.], 
        rows: 3, 
        cols: 3
    };
    println!("A = {:?}", square);
    let (plu, permut) = square.decompose_lup(0.0001).unwrap();
    assert_eq!(plu, 
        Matrix { elem: vec![10.0, -9.0, 18.0, 
            0.7, -5.7000003, -1.5999994, 
            -1.0, -0.17543858, 20.7193], 
            rows: 3, cols: 3 });
    println!("LU = {:?} \n P = {:?}", plu, permut);

    let inverse = plu.inverse(&permut);
    match inverse{
        Ok(inverse) => {
            println!("inv A = {:?}", inverse);
            assert_eq!(inverse, 
                Matrix { elem: 
                    vec![0.12362405, -0.17527518, -0.09906857, 
                    0.11092294, -0.1778154, -0.013547834, 
                    0.042337, 0.0084673995, 0.04826418], 
                    rows: 3, 
                    cols: 3 })
        },
        Err(e) => panic!(e)
    }

    let det = plu.determinant(&permut);
    assert_eq!(det, -1181.0001);
    println!("det A = {:?}", det);

    let b = Matrix{
        elem: vec![1.,2.,3.], 
        rows: 3, 
        cols: 1
    };

    println!("Solving AX = [1,2,3]");
    let x = plu.solve_lup(&permut, &b).unwrap();
    assert_eq!(x, Matrix{elem: vec![-0.524132, -0.28535134, 0.20406434], 
        rows: 3, cols: 1});
    println!("X = {:?}", x.elem);
    let verif = square.dot(&x);
    match verif{
        Ok(verif) => {
            assert_eq!(verif, Matrix{elem: vec![1.0, 1.9999998, 3.0], 
                rows: 3, cols: 1});
            println!("Verification: A.x = {:?}", square.dot(&x));
        },
        Err(e) => panic!(e)
    }

}

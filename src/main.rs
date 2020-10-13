mod rmatrix;

fn main() {

    let mut s1 = rmatrix::RMatrix{elem: vec![1.,-2.,2.,-1.], rows: 2, cols: 2};
    let s2 = rmatrix::RMatrix{elem: vec![3.,4.,-1.,2.], rows: 2, cols: 2};
    println!("s1 = {:?}", s1);
    println!("s2 = {:?}", s2);

    let s3 = s1.add(&s2);
    match s3{
        Ok(s3) => println!("s1+s2 = {:?}", s3),
        Err(e) => panic!(e)
    }

    s1.mut_add(&s2);
    println!("s1+s2 = {:?}", s1);

    let s4 = rmatrix::RMatrix{elem: vec![1.,3.,5.,2.,4.,6.], rows: 2, cols: 3};
    let s5 = rmatrix::RMatrix{elem: vec![-1.,1.,2.,3.,0.,-2.], rows: 3, cols: 2};
    println!("s4 = {:?}", s4);
    println!("s5 = {:?}", s5);

    let mut s6 = s4.dot(&s5);
    match s6{
        Ok(s6) => println!("s4xs5 = {:?}", s6),
        Err(e) => panic!(e)
    }
    s6 = s5.dot(&s4);
    match s6{
        Ok(s6) => println!("s5xs4 = {:?}", s6),
        Err(e) => panic!(e)
    }

    println!("s4' = {:?}", s4.transpose());
    println!("s5' = {:?}", s5.transpose());

    // LU decomposition, inversion, solving, determinant
    let mut square = rmatrix::RMatrix{
        elem: vec![10.,-9.,18.,7.,-12., 11.,-10.,10.,3.], 
        rows: 3, 
        cols: 3
    };
    println!("A = {:?}", square);
    let (plu, permut) = square.decompose_lup(0.0001).unwrap();
    println!("LU = {:?} \n P = {:?}", plu, permut);

    let inverse = plu.inverse(&permut);
    match inverse{
        Ok(inverse) => println!("inv A = {:?}", inverse),
        Err(e) => panic!(e)
    }

    let det = plu.determinant(&permut);
    println!("det A = {:?}", det);

    let b = rmatrix::RMatrix{
        elem: vec![1.,2.,3.], 
        rows: 3, 
        cols: 1
    };

    println!("Solving AX = [1,2,3]");
    let x = plu.solve_lup(&permut, &b).unwrap();
    println!("X = {:?}", x.elem);
    println!("Verification: A.x = {:?}", square.dot(&x));
}


// 定义函数模板
template<int M, int N>
LinearProgrammingResult<M, N> BoundedRevisedSimplex(LinearProgrammingProblem<M, N> problem) {
    LinearProgrammingResult<M, N> result;
    // 实现线性规划算法
    // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e
    float tol=1e-7;
    const int n_m=N-M;
    int* nind = generateSequence(0, n_m-1);
    std::cout << "nind: [";
    for (size_t i = 0; i <N- M; ++i) {
        std::cout << nind[i];
        if (i <N- M - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    int* ind_all = generateSequence(0, N-1);
    std::cout << "ind_all: [";
    for (size_t i = 0; i <N; ++i) {
        std::cout << ind_all[i];
        if (i <N- 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    setdiff(ind_all, N, problem.inB, M, problem.inD);
    std::cout << "problem.inB: [";
    for (size_t i = 0; i <M; ++i) {
        std::cout << problem.inB[i];
        if (i <M- 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << "problem.inD: [";
    for (size_t i = 0; i <N-M; ++i) {
        std::cout << problem.inD[i];
        if (i <N-M- 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;


    std::cout << "A:" << std::endl;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            std::cout << problem.A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "b: [";
    for (size_t i = 0; i < M; ++i) {
        std::cout << problem.b[i];
        if (i < M - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "c: [";
    for (size_t i = 0; i < N; ++i) {
        std::cout << problem.c[i];
        if (i < N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "h: [";
    for (size_t i = 0; i < N; ++i) {
        std::cout << problem.h[i];
        if (i < N- 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    
    //

    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if(!problem.e[j])
            {
                problem.A[i][j] *=-1;
                problem.b[i]+=problem.A[i][j]*problem.h[j];
            }
        }
    }
    for(int j=0; j<N; ++j)
    {
        if(!problem.e[j])
        {
            problem.c[j] *=-1;
        }
    }
    std::cout << "A:" << std::endl;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            std::cout << problem.A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "b: [";
    for (size_t i = 0; i < M; ++i) {
        std::cout << problem.b[i];
        if (i < M - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "c: [";
    for (size_t i = 0; i < N; ++i) {
        std::cout << problem.c[i];
        if (i < N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "h: [";
    for (size_t i = 0; i < N; ++i) {
        std::cout << problem.h[i];
        if (i < N- 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    
    //==============================
    matrix::SquareMatrix<float, M> A_inB;
    matrix::Matrix<float, M, n_m> A_inD;
    matrix::Vector<float, M> c_inB;
    matrix::Vector<float, n_m> c_inD;
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<M; ++j)
        {
            A_inB(i,j)=problem.A[i][problem.inB[j]];
            if(j<n_m)
            {
                A_inD(i,j)=problem.A[i][problem.inD[j]];
            }
        }
        c_inB(i)=problem.c[problem.inB[i]];
    }
    for(int i=0; i<n_m; ++i)
    {
        c_inD(i)=problem.c[problem.inD[i]];
    }
    matrix::Vector<float, M> b_vec(problem.b);

    // inital some value
    Matrix<float, 1, M> lamt;
    lamt.setZero();
    Matrix<float, 1, n_m> rdt;
    rdt.setZero();
    matrix::Vector<float, M> A_qel;
    A_qel.setZero();
    matrix::Vector<float, M> yq;
    yq.setZero();
    matrix::Vector<float, M> rat;
    rat.setZero();
    
    //  %Initial Solution
    matrix::Vector<float, M> y0 = inv(A_inB)*b_vec;
    std::cout << "y0:";
    y0.print();
    bool done = false;
    bool unbounded = false;
    int iters =0;
     while ((!done  || !unbounded ) && (iters <= problem.itlim))
    {
        iters = iters+1;
        // lamt= (inv(A_inB).transpose()*c_inB).transpose();
        matrix::LeastSquaresSolver<float, M,M> LSsolver_lamt(A_inB.transpose());
        lamt = LSsolver_lamt.solve(c_inB).transpose();
        rdt = c_inD.transpose()-lamt*A_inD;
        std::cout << "lamt:";
        lamt.print();
        std::cout << "rdt:";
        rdt.print();
        float minr;
        size_t qind;
        min(rdt.transpose(), minr, qind);
        std::cout << "minr:"<<minr<<std::endl;
        std::cout << "qind:"<<qind<<std::endl;
        if(minr >=0)  // If all relative costs are positive then the solution is optimal
        { 
            done = true;
            break;
        }
        int qel = problem.inD[qind];  // Unknown to Enter the basis minimizes relative cost
        A_qel(0)=problem.A[0][qel];
        A_qel(1)=problem.A[1][qel];
        A_qel(2)=problem.A[2][qel];
        std::cout << "qel:"<<qel<<std::endl;

        // yq=inv(A_inB)* A_qel; // Vector to enter in terms of the current Basis vector
        matrix::LeastSquaresSolver<float, M,M> LSsolver1(A_inB);
        yq = LSsolver1.solve(A_qel);
        std::cout << "yq:";
        yq.print();
        bool flag=false;
        
        for(int i=0;i<M;++i){
            if(std::abs(yq(i)) > tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag)
        {
            unbounded = true; // Check this condition
            break;
        }
        // Recompute rations and determine variable to leave
        
        float hinB[M];
        for(int i=0;i<M;++i)
        {
            if(std::abs(yq(i))>tol)
            {
                rat(i)=y0(i)/yq(i);
                
            }
            else
            {
                rat(i)=INFINITY;
                /* code */
            }
        }
        for(int i=0;i<M;++i)
        {
            hinB[i]=problem.h[problem.inB[i]];
            if(yq(i)<0 && std::abs(yq(i))>tol )
            {                    
                rat(i)-=hinB[i]/yq(i);
            }
        }
        std::cout << "rat:";
        rat.print();
         // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        float minrat=rat(0);
        size_t p=0;
        min(rat, minrat, p);
        std::cout << "minrat:"<<minrat<<std::endl;
        std::cout << "p:"<<p<<std::endl;
        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if (std::abs(minrat) <= tol)
        {
            //Find negative relative cost
            std::cout << " Find negative relative cost "<< std::endl; 
            for(int i=0;i<N-M;++i)
            {
                // std::cout << "rdt(0,i):"<<rdt(0,i)<<std::endl; 
                if(rdt(0,i)<0){ //Note that since minr <0 indm is not empty 
                    qind=nind[i];
                    qel = problem.inD[qind];//Unknown to Enter the basis is first indexed to avoid cycling
                    std::cout << "qind:"<<qind<<std::endl;
                    std::cout << "qel:"<<qel<<std::endl;
                    break;
                }
            }
            A_qel(0)=problem.A[0][qel];
            A_qel(1)=problem.A[1][qel];
            A_qel(2)=problem.A[2][qel];
            // yq=inv(A_inB)* A_qel;
            matrix::LeastSquaresSolver<float, M,M> LSsolver2(A_inB);
            yq = LSsolver2.solve(A_qel);
            std::cout << "yq:";
            yq.print();
            bool flag=false;
            for(int i=0;i<M;++i){
                if(std::abs(yq(i)) > tol)
                {
                    flag = true; // Check this condition
                    break;
                }
            }
            if(!flag)
            {
                unbounded = true; // Check this condition
                // break;
            }
            // Recompute rations and determine variable to leave
            // Recompute rations and determine variable to leave
            float hinB[M];
            for(int i=0;i<M;++i)
            {
                hinB[i]=problem.h[problem.inB[i]];
                if(std::abs(yq(i))>tol)
                {
                    rat(i)=y0(i)/yq(i);
                    if(yq(i)<0)
                    {                    
                        rat(i)-=hinB[i]/yq(i);
                    }
                }
                else
                {
                    rat(i)=INFINITY;
                    /* code */
                }
            }
            std::cout << "rat:";
            rat.print();
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            minrat=rat(0);
            p=0;
            min(rat, minrat, p);
            std::cout << "minrat:"<<minrat<<std::endl;
            std::cout << "p:"<<p<<std::endl;
        }
        if (minrat >= problem.h[qel])
        {
            std::cout << " Case 1 "<< std::endl; 
            problem.e[qel] =!problem.e[qel];
            for(int i=0; i<M; ++i)
            {
                problem.A[i][qel] *= -1;
                b_vec(i)+=problem.A[i][qel]*problem.h[qel];
            }
            problem.c[qel] *= -1;

        }
        else if(yq(p) > 0)
        {
            std::cout << " Case 21 "<< std::endl; 
            int pel = problem.inB[p];
            std::cout << "pel:"<<pel<<std::endl;
            std::cout << "qel:"<<qel<<std::endl;
            problem.inB[p]= qel;
            problem.inD[qind]= pel;
            // update x_inX
            for(int i=0; i<M; ++i)
            {
                A_inB(i,p)=problem.A[i][qel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inB(p)=problem.c[qel];
            }
            for(int i=0; i<M; ++i)
            {
                A_inD(i,qind)=problem.A[i][pel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inD(qind)=problem.c[pel];
            }
            std::cout << "A:" << std::endl;
            for (size_t i = 0; i < M; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    std::cout << problem.A[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "b: [";
            for (size_t i = 0; i < M; ++i) {
                std::cout << problem.b[i];
                if (i < M - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
            std::cout << "b_vec:";
            b_vec.print();
            std::cout << "problem.inB: [";
            for (size_t i = 0; i <M; ++i) {
                std::cout << problem.inB[i];
                if (i <M- 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
            std::cout << "problem.inD: [";
            for (size_t i = 0; i <N-M; ++i) {
                std::cout << problem.inD[i];
                if (i <N-M- 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }
        else
        {
            std::cout << " Case 22 "<< std::endl; 
            int pel = problem.inB[p];
            std::cout << "pel:"<<pel<<std::endl;
            problem.e[pel]=!problem.e[pel];
            std::cout << "problem.e[pel]:"<<problem.e[pel]<<std::endl;
            for(int i=0; i<M; ++i)
            {
                problem.A[i][pel] *= -1;
                b_vec(i)+=problem.A[i][pel]*problem.h[pel];
                problem.b[i]=b_vec(i);

            }
            problem.inB[p]= qel;
            problem.inD[qind]= pel;
            problem.c[pel] *= -1;
            // update x_inX
            for(int i=0; i<M; ++i)
            {
                A_inB(i,p)=problem.A[i][qel];
            }
            
            for(int i=0; i<n_m; ++i)
            {
                c_inB(p)=problem.c[qel];
            }
            for(int i=0; i<M; ++i)
            {
                A_inD(i,qind)=problem.A[i][pel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inD(qind)=problem.c[pel];
            }
            std::cout << "A:" << std::endl;
            for (size_t i = 0; i < M; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    std::cout << problem.A[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "b: [";
            for (size_t i = 0; i < M; ++i) {
                std::cout << problem.b[i];
                if (i < M - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
            std::cout << "b_vec:";
            b_vec.print();
            std::cout << "c: [";
            for (size_t i = 0; i < N; ++i) {
                std::cout << problem.c[i];
                if (i < N - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;

            std::cout << "problem.inB: [";
            for (size_t i = 0; i <M; ++i) {
                std::cout << problem.inB[i];
                if (i <M- 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
            std::cout << "problem.inD: [";
            for (size_t i = 0; i <N-M; ++i) {
                std::cout << problem.inD[i];
                if (i <N-M- 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }
        // y0=inv(A_inB)* b_vec;

        matrix::LeastSquaresSolver<float, M,M> LSsolver(A_inB);
        y0 = LSsolver.solve(b_vec);
        std::cout << "y0:";
        y0.print();
        std::cout << "A_inB:";
        A_inB.print();
        





        
    }
    result.errout = unbounded; 
    // y0.print();
    // 设置 result.y0, result.inB, result.e 等结果
    for(int i=0; i<M; ++i)
    {
        result.y0[i]=y0(i);
        result.inB[i]=problem.inB[i];
    }
    for(int i=0; i<N; ++i)
    {
        result.e[i]=problem.e[i];
    }
    result.iters=iters;
    std::cout << "iters:"<<iters<<std::endl;
    std::cout << "problem.itlim:"<<problem.itlim<<std::endl;
    return result;
}


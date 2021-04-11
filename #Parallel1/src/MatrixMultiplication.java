import java.time.Duration;
import java.time.Instant;
import java.util.Scanner;

public class MatrixMultiplication {
	public static int M;
    public static int K;
    public static int N;

    public static void main(String[] args) {  
    	Instant start = Instant.now();
    	Scanner input = new Scanner(System.in);
		System.out.println("Enter the number M: ");
		M = input.nextInt();
		System.out.println("Enter the number K: ");
		K = input.nextInt();
		System.out.println("Enter the number N: ");
		N = input.nextInt();
		
		WorkerThread [][] Threads = new WorkerThread [M][N];
		 
		int[][] A = new int[M][K];
		int[][] B = new int[K][N];
		int[][] C = new int[M][N];
		
		System.out.println("Enter the numbers of the first matrix.\n");
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < K; j++) {
				A[i][j] = input.nextInt();
			}
		}
		System.out.println("Enter the numbers of the 2nd matrix.\n");
		for (int i = 0; i <K; i++) {
			for (int j = 0; j < N; j++) {
				B[i][j] = input.nextInt();
			}
		}
		
        //creates 9 Worker threads. Each thread Calculates a Matrix Value and sets it to C matrix
        for (int i = 0; i<M; i++){
            for (int j=0; j<N; j++){
                Threads[i][j] = new WorkerThread(i,j,A,B,C);
                Threads[i][j].run();
                Threads[i][j].start();
            }
        }

        //Outputs the Values of Matrix C
        System.out.println("Elements of Matrix C:");
        for (int i = 0; i<M; i++){
            for (int j=0; j<N; j++){
                System.out.println("["+i+","+j+"] = "+ C[i][j]);
            }
        }  
        

		Instant end = Instant.now();

		Duration interval = Duration.between(start, end);

		System.out.println("Execution time in seconds: " + interval.getSeconds());
    }
} 
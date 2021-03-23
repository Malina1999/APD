import java.time.Duration;
import java.time.Instant;
import java.util.Scanner;

public class Multiplication {
	public static void main(String args[]) throws InterruptedException {
		Instant start = Instant.now();
		int n;
		Scanner input = new Scanner(System.in);
		System.out.println("Enter the number of rows and columns of the matrices. They must be equal.");
		n = input.nextInt();
		int[][] a = new int[n][n];
		int[][] b = new int[n][n];
		int[][] c = new int[n][n];
		System.out.println("Enter the numbers of the first matrix.\n");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				a[i][j] = input.nextInt();
			}
		}
		System.out.println("Enter the numbers of the 2nd matrix.\n");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				b[i][j] = input.nextInt();
			}
		}
		System.out.println("Generating the multiplication of matrices.....");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					c[i][j] = c[i][j] + a[i][k] * b[k][j];
				}
			}
		}
		System.out.println("The product of the matrices is shown as below");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				System.out.print(c[i][j] + " ");
			}
			System.out.println();
		}
		input.close();

		Instant end = Instant.now();

		Duration interval = Duration.between(start, end);

		System.out.println("Execution time in seconds: " + interval.getSeconds());

	}
}
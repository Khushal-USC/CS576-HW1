
import java.awt.*;
import java.awt.image.*;
import java.io.*;

import javax.imageio.ImageIO;
import javax.swing.*;


public class ImageDisplay {

	JFrame frame;
	JLabel lbIm1;
	BufferedImage imgOne;

	// Modify the height and width values here to read and display an image with
  	// different dimensions. 
	

	// int width = 256;
	// int height = 170;

	// int width = 128;
	// int height = 128;
	

	/** Read Image RGB
	 *  Reads the image of given width and height at the given imgPath into the provided BufferedImage.
	 */
	private void readImageRGB(int width, int height, String imgPath, BufferedImage img)
	{
		try
		{
			int frameLength = width*height*3;

			File file = new File(imgPath);
			RandomAccessFile raf = new RandomAccessFile(file, "r");
			raf.seek(0);

			long len = frameLength;
			byte[] bytes = new byte[(int) len];

			raf.read(bytes);

			int ind = 0;
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width; x++)
				{
					byte a = 0;
					
					byte r = bytes[ind];
					byte g = bytes[ind+(height*width)];
					byte b = bytes[ind+(height*width*2)]; 

					int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
					//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
					//int pix = 0xff32A852;
					img.setRGB(x,y,pix);
					ind++;
				}
			}
		}
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	private void generateKernal(double rho , double[][] matrix)
	{
		int midX = matrix.length/2, midY = matrix.length/2;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				int x = i - midX;
				int y = j - midY;
				double g = (1/(2*Math.PI*Math.pow(rho, 2))) * Math.pow(Math.E, -(((x*x) + (y*y))/(2*Math.pow(rho, 2))));
				matrix[i][j] = g;
			}
		}
	}

	private void gaussBlur(int width, int height, int kernalSize, double rho, BufferedImage imgIn, BufferedImage imgOut)
	{

		double[][] m = new double[kernalSize][kernalSize];
		generateKernal(rho,m);
		int R = kernalSize/2;
		//System.out.println("AASDASD" + R);
		// for (double[] ds : m) {
		// 	for (double ds2 : ds) {
		// 		System.out.print(ds2 + " ");
		// 	}
		// 	System.out.println();
		// }

		for (int i = R; i < height-R; i++) {
			// System.err.println("herey" + i);
			for (int j = R; j < width-R; j++) {
				// System.err.println("here" + j);
				double r_s = 0;
				double g_s = 0;
				double b_s = 0;
				double w_s = 0;
				//byte a = (byte)((color & 0xff000000) >> 24);
				int y = 0;
				for (int k = i-R; k <= i+R; k++) {
					int x = 0;
					for (int l = j-R; l <= j+R; l++) {
						int color = imgIn.getRGB(l, k);
						double weight = m[x][y];
						// System.err.println(weight);

						// Components will be in the range of 0..255:
						r_s += (((color & 0xff0000) >> 16) * weight);
						g_s += (((color & 0xff00) >> 8) * weight);
						b_s += ((color & 0xff) * weight);
						w_s += weight;
						
						x++;
						
					}
					y++;
				}
				// System.out.println(r_s + " " + b_s + " " + g_s + " ");
				// System.out.println(w_s);
				//System.out.println(count);
				
				byte r = (byte)(r_s / w_s);
				byte b = (byte)(b_s / w_s);
				byte g = (byte)(g_s / w_s);
				// r = (byte)Math.min(255, Math.max(0, r));
				// g = (byte)Math.min(255, Math.max(0, g));
				// b = (byte)Math.min(255, Math.max(0, b));
				//System.out.println(r + " " + b + " " + g + " ");

				int pixOut = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				imgOut.setRGB(j,i,pixOut);
			}
		}

	}


	private void downSample1(int width, int height,int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut)
	{
		int y = 0;
		int outY = 0;
		while(outY < outHeight)
		{
			int x = 0;
			int outX = 0;
			while(outX < outWidth)
			{
				//System.out.println("herex " +x);
				int color = imgIn.getRGB(x, y);

				// Components will be in the range of 0..255:
				byte r = (byte)((color & 0xff0000) >> 16);
				byte g = (byte)((color & 0xff00) >> 8);
				byte b = (byte)(color & 0xff);
				
				byte a = (byte)((color & 0xff000000) >> 24);

				//int pixOut = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				int pix = ((a << 24) + (r << 16) + (g << 8) + b);
				//int pix = 0xff32A852;
				imgOut.setRGB(outX,outY,pix);
				x+=(width/outWidth);
				outX++;
			}
			//System.out.println("herey " +y);
			y+=(height/outHeight);
			outY++;
		}
	}

	static BufferedImage deepCopy(BufferedImage bi) {
		ColorModel cm = bi.getColorModel();
		boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
		WritableRaster raster = bi.copyData(null);
		return new BufferedImage(cm, raster, isAlphaPremultiplied, null);
	   }


	public void showIms(String[] args){
		int width = 400;
		int height = 300;

		// Read in the specified image
		imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		readImageRGB(width, height, args[0], imgOne);

		// //Downsample1
		// int outWidth = 640;
		// int outHeight = 400;
		// BufferedImage imgDS1Out = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		// downSample1(width, height, outWidth, outHeight, imgOne, imgDS1Out);

		//GaussBlur
		double rho = 3;
		
		BufferedImage imgBlur = deepCopy(imgOne);
		gaussBlur(width, height, (int)(3*rho), rho, imgOne, imgBlur);

		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		// lbIm1 = new JLabel(new ImageIcon(imgOne));
		// lbIm1 = new JLabel(new ImageIcon(imgDS1Out));
		lbIm1 = new JLabel(new ImageIcon(imgBlur));

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		frame.getContentPane().add(lbIm1, c);

		frame.pack();
		frame.setVisible(true);

		File outputfile = new File("output2.png");
		// try {
		// 	ImageIO.write(imgOne, "png", outputfile);
		// } catch (IOException e) {
		// 	// TODO Auto-generated catch block
		// 	e.printStackTrace();
		// }
		// try {
		// 	ImageIO.write(imgDS1Out, "png", outputfile);
		// } catch (IOException e) {
		// 	// TODO Auto-generated catch block
		// 	e.printStackTrace();
		// }
		try {
			ImageIO.write(imgBlur, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		ImageDisplay ren = new ImageDisplay();
		ren.showIms(args);
		// double[][] m = new double[3][3];
		System.out.println();
		// ren.generateKernal(3,1,m);
		// for (double[] ds : m) {
		// 	for (double ds2 : ds) {
		// 		System.out.print(ds2 + " ");
		// 	}
		// 	System.out.println();
		// }
	}

}

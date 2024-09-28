
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
		// System.out.println("AASDASD" + R);
		// for (double[] ds : m) {
		// 	for (double ds2 : ds) {
		// 		System.out.print(ds2 + " ");
		// 	}
		// 	System.out.println();
		// }

		for (int i = 0; i < height; i++) {
			// System.err.println("herey" + i);
			for (int j = 0; j < width; j++) {
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
						int color = 0;
						if(l >= 0 && l < width && k >=0 && k < height){
							color = imgIn.getRGB(l, k);
						} else {
							color = imgIn.getRGB(Math.min(Math.max(l,0),width-1), Math.min(Math.max(k,0),height-1));
						}
						
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

	private int linInterp(double srcX, double srcY, int x1, int y1, int x2, int y2, int q11, int q12, int q21, int q22){
		//X dir interp
		int r1_x_y = 0;
		if (x1 == x2) {
			r1_x_y = q11;  // No interpolation in x direction, use q11 directly
		} else {
			r1_x_y = (int)(q11 *((double)(x2 - srcX)/(x2-x1)) + q21*((double)(srcX-x1)/(x2-x1)));
		}

		double r2_x_y;
		if (x1 == x2) {
			r2_x_y = q12;  // No interpolation in x direction, use q12 directly
		} else {
			r2_x_y = (int)(q12 *((double)(x2 - srcX)/(x2-x1)) + q22*((double)(srcX-x1)/(x2-x1)));
		}

		double p_x_y;
		if (y1 == y2) {
			p_x_y = r1_x_y;  // No interpolation in y direction, use r1 directly
		} else {
			p_x_y = (int)(r1_x_y *((double)(y2 - srcY)/(y2-y1)) + r2_x_y*((double)(srcY-y1)/(y2-y1)));
		}
		return (int)Math.min(Math.max(p_x_y,0),255);
	}

	private int getR(int color){
		return ((color & 0xff0000) >> 16);
	}
	private int getG(int color){
		return ((color & 0xff00) >> 8);
	}
	private int getB(int color){
		return (color & 0xff);
	}

	private void bilinear(int width, int height,int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut)
	{
		System.out.println("doing 2");
		int outY = 0;
		while(outY < outHeight){
			int outX = 0;
			while(outX < outWidth)
			{
				double scaleX = (double) width / outWidth;
				double scaleY = (double) height / outHeight;

				double srcX = (outX * scaleX);
				double srcY = (outY * scaleY);
				int x1 = (int)srcX;
				int x2 = Math.min(Math.max(x1 + 1,0),width-1);
				int y1 = (int)srcY;
				int y2 = Math.min(Math.max(y1 + 1,0),height-1);

				int q11 = imgIn.getRGB(x1, y1);
				int q12 = imgIn.getRGB(x1, y2);
				int q21 = imgIn.getRGB(x2, y1);
				int q22 = imgIn.getRGB(x2, y2);
				

				// Components will be in the range of 0..255:
				byte r = (byte)linInterp(srcX, srcY, x1, y1, x2, y2, getR(q11), getR(q12), getR(q21), getR(q22));
				byte g = (byte)linInterp(srcX, srcY, x1, y1, x2, y2, getG(q11), getG(q12), getG(q21), getG(q22));
				byte b = (byte)linInterp(srcX, srcY, x1, y1, x2, y2, getB(q11), getB(q12), getB(q21), getB(q22));
				
				//byte a = (byte)((color & 0xff000000) >> 24);

				int pixOut = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
				//int pix = 0xff32A852;
				imgOut.setRGB(outX,outY,pixOut);
				outX++;
			}
			//System.out.println("herey " +y);
			outY++;
		}
	}

	private void upSample1(int width, int height,int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut)
	{
		System.out.println("doing 1");
		int outY = 0;
		while(outY < outHeight)
		{
			int outX = 0;
			while(outX < outWidth)
			{
				float scaleX = (float) width / outWidth;
				float scaleY = (float) height / outHeight;

				int srcX = (int) (outX * scaleX);
				int srcY = (int) (outY * scaleY);
				int color = imgIn.getRGB(Math.min(width - 1, srcX), Math.min(height - 1, srcY));

				// Components will be in the range of 0..255:
				byte r = (byte)((color & 0xff0000) >> 16);
				byte g = (byte)((color & 0xff00) >> 8);
				byte b = (byte)(color & 0xff);
				
				byte a = (byte)((color & 0xff000000) >> 24);

				//int pixOut = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				int pix = ((a << 24) + (r << 16) + (g << 8) + b);
				//int pix = 0xff32A852;
				imgOut.setRGB(outX,outY,pix);
				outX++;
			}
			//System.out.println("herey " +y);
			outY++;
		}
	}


	private void downSample1(int width, int height,int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut)
	{
		int y = 0;
		int outY = 0;
		while(outY < outHeight/3)
		{
			int x = 0;
			int outX = 0;
			while(outX < outWidth/4)
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

	// private double nonLinearStep(int width, int outWidth, int step){
	// 	return (width * Math.pow(step, 2))/(Math.pow(outWidth, 2));
	// }
	// private int linearStep(int width, int outWidth, int step){
	// 	return (outWidth/width) * step;
	// }

	private int linearStep(int width, int outWidth, int step) {
		return (int)((double)step * width / outWidth);
	}
	private int linearStep2(int width, int outWidth, int step) {
		return (int)(((double)step * width / outWidth)-((0.0004)*step*step));
	}

	private int nonLinearStep(int width, int outWidth, int step){
		return (int)((width * Math.pow(step, 2))/(Math.pow(outWidth, 2)));
	}
	private int nonLinearStepFixed(int width,int height, int outWidth, int step){
		// return (int)((width * Math.pow(step, 0.6))/((Math.pow(outWidth, 0.6))));
		return (int)(((double)step * width / outWidth)+((0.003)*step*step));
	}

	private int nonLinearStepFixed2(int width, int height, int outWidth, int step){
		// return (int)((width * Math.pow(step, 0.6))/((Math.pow(outWidth, 0.6))));
		return (int)(((double)step * width*5 / height)+((0.0003)*step*step));
	}
	// private int nonLinearStepFixed(int width, int outWidth, int step) {
	// 	// Reverse the stretching effect: distance from the edge instead of the center
	// 	int maxStep = outWidth - 1;  // Maximum value for step
	// 	int invertedStep = maxStep - step;  // Invert the step to stretch more at the edges
	
	// 	// Apply the cubic scaling to the inverted step
	// 	return (int)((width * Math.pow(invertedStep, 3)) / (2 * Math.pow(maxStep, 3)));
	// }
	

	public static int clamp(int value, int min, int max) {
		return (int)Math.max(min, Math.min(max, value));
	}

	// public static BufferedImage scaleImage(BufferedImage originalImage, int newWidth, int newHeight) {
    //     // Create a new buffered image with the desired width and height
    //     BufferedImage scaledImage = new BufferedImage(newWidth, newHeight, BufferedImage.TYPE_INT_ARGB);
        
    //     // Get the Graphics2D object from the new image
    //     Graphics2D g2d = scaledImage.createGraphics();
        
    //     // Set rendering hints for better quality
    //     g2d.setRenderingHint(java.awt.RenderingHints.KEY_INTERPOLATION, java.awt.RenderingHints.VALUE_INTERPOLATION_BILINEAR);
    //     g2d.setRenderingHint(java.awt.RenderingHints.KEY_RENDERING, java.awt.RenderingHints.VALUE_RENDER_QUALITY);
    //     g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
        
    //     // Draw the original image scaled to the new dimensions
    //     g2d.drawImage(originalImage, 0, 0, newWidth, newHeight, null);
        
    //     // Dispose the graphics context to free up resources
    //     g2d.dispose();
        
    //     return scaledImage;
    // }

	public static BufferedImage nonlinearHorizontalScale(BufferedImage originalImage, int newWidth, double stretchFactor) {
        int originalWidth = originalImage.getWidth();
        int originalHeight = originalImage.getHeight();

        // Create a new buffered image with the desired width and the same height
        BufferedImage scaledImage = new BufferedImage(newWidth, originalHeight, BufferedImage.TYPE_INT_ARGB);

        // Loop over the new image pixels
        for (int y = 0; y < originalHeight; y++) {
            for (int x = 0; x < newWidth; x++) {
                // Calculate the corresponding original x position based on a nonlinear transformation
                double normalizedX = (double) x / newWidth; // Normalize x to a value between 0 and 1
                double transformedX = nonlinearTransformation(normalizedX, stretchFactor);

                int originalX = (int) (transformedX * originalWidth); // Get the corresponding original x position

                // Bound the originalX to be within the original image's dimensions
                originalX = Math.max(0, Math.min(originalWidth - 1, originalX));

                // Get the color of the original pixel
                Color color = new Color(originalImage.getRGB(originalX, y));

                // Set the color in the new image
                scaledImage.setRGB(x, y, color.getRGB());
            }
        }

        return scaledImage;
    }

    // Nonlinear transformation function
    private static double nonlinearTransformation(double x, double stretchFactor) {
        // Stretch more near the left or right edge based on the stretch factor
        if (stretchFactor > 0) {
            // Stretch more near the right edge
            return Math.pow(x, 1 / (1 + stretchFactor));
        } else {
            // Stretch more near the left edge
            return 1 - Math.pow(1 - x, 1 / (1 - stretchFactor));
        }
    }


	private void downSampleFirst(int width, int height, int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut) {
		// Calculate aspect ratios
		double inputAspect = (double) width / height;
		double outputAspect = (double) outWidth / outHeight;
		
		int newWidth, newHeight;
		
		// If the input image is wider than the output, scale based on width

		// stretchVert is true if the input image is taller than the output and letterboxing will occur on the top and bottom
		boolean stretchVert;
		if (inputAspect > outputAspect) {
			stretchVert = true;
			newWidth = outWidth;
			newHeight = (int) (outWidth / inputAspect);
		} else {
			stretchVert = false;
			// Otherwise, scale based on height
			newHeight = outHeight;
			newWidth = (int) (outHeight * inputAspect);
		}
		
		// Create a scaled version of the input image
		BufferedImage scaledImage = new BufferedImage(newWidth, newHeight, imgIn.getType());
		Graphics2D g2d = scaledImage.createGraphics();
		g2d.drawImage(imgIn, 0, 0, newWidth, newHeight, null);
		g2d.dispose();
		
		// Fill the output image with black (or any background color) to handle letterboxing
		Graphics2D g2 = imgOut.createGraphics();
		g2.setColor(Color.BLACK); // Set letterbox color
		g2.fillRect(0, 0, outWidth, outHeight); // Fill the entire output image
		
		// Calculate position to center the scaled image
		int xOffset = (outWidth - newWidth) / 2;
		int yOffset = (outHeight - newHeight) / 2;

		//Get coordinates of center of imgOut
		int centerX = outWidth/2;
		int centerY = outHeight/2;

		System.out.println();
		System.out.println("centerX: " + centerX);
		System.out.println("centerY: " + centerY);
		//print newWidth and newHeight
		System.out.println("newWidth: " + newWidth);
		System.out.println("newHeight: " + newHeight);

		// Stretch pixels outward from edges to fill letterbox areas

		// R represents the distance from the center of the image whic
		int left_edgeX = centerX - newWidth/2;
		int right_edgeX = centerX + newWidth/2;

		//pront left_edgeX and right_edgeX
		System.out.println("left_edgeX: " + left_edgeX);
		System.out.println("right_edgeX: " + right_edgeX);

		int r = 320;

		//Get difference between newWidth and outWidth
		int diffWidth = outWidth - newWidth;

		int leftSampleX = left_edgeX + r;
		int rightSampleX = right_edgeX - r;

		//Get the subimage of the scaled image using edgeX and r
		BufferedImage subImageLeft = scaledImage.getSubimage(0, 0, r, outHeight);
		//Scale the subimageleft to have width of diffWidth
		subImageLeft = nonlinearHorizontalScale(subImageLeft, leftSampleX, 0);

		//Draw the scaled subimage to the left of the output image
		//g2.drawImage(subImageLeft, 0, 0, null);

		//Get the subimage of the scaled image using edgeX and r
		BufferedImage subImageRight = scaledImage.getSubimage(newWidth-r, 0, r, outHeight);
		//Scale the subimage to have width of diffWidth
		subImageRight= nonlinearHorizontalScale(subImageRight, rightSampleX, 0);

		// //Output the subimage to the right of the output image
		// g2.drawImage(subImageRight, outWidth-r, 0, null);

		//Draw the scaled subimage to the right of the output image
		

		// Draw the scaled image onto the output image
		g2.drawImage(scaledImage, xOffset, yOffset, null);
		g2.drawImage(subImageLeft, 0, 0, null);
		g2.drawImage(subImageRight, rightSampleX, 0, null);

		//Output the subimageright to file right.png
		File outputfile = new File("right.png");
		try {
			ImageIO.write(subImageRight, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		g2.dispose();
	}

	private void downSampleFUCK(int width, int height, int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut, double stretchRadius) {
		// Calculate aspect ratios
		double inputAspect = (double) width / height;
		double outputAspect = (double) outWidth / outHeight;
		
		int newWidth, newHeight;
		
		// If the input image is wider than the output, scale based on width
		if (inputAspect > outputAspect) {
			newWidth = outWidth;
			newHeight = (int) (outWidth / inputAspect);
		} else {
			// Otherwise, scale based on height
			newHeight = outHeight;
			newWidth = (int) (outHeight * inputAspect);
		}
		
		// Create a scaled version of the input image
		BufferedImage scaledImage = new BufferedImage(newWidth, newHeight, imgIn.getType());
		Graphics2D g2d = scaledImage.createGraphics();
		g2d.drawImage(imgIn, 0, 0, newWidth, newHeight, null);
		g2d.dispose();
		
		// Fill the output image with the original image
		Graphics2D g2 = imgOut.createGraphics();
		
		// Calculate position to center the scaled image
		int xOffset = (outWidth - newWidth) / 2;
		int yOffset = (outHeight - newHeight) / 2;
		
		// Copy the scaled image to the center of the output
		g2.drawImage(scaledImage, xOffset, yOffset, null);
		
		// Stretch pixels outward from edges to fill letterbox areas
		for (int y = 0; y < outHeight; y++) {
			for (int x = 0; x < outWidth; x++) {
				// Check if we're in the letterbox areas
				if (x < xOffset || x >= xOffset + newWidth || y < yOffset || y >= yOffset + newHeight) {
					// Calculate distance from center of the image
					double distX = Math.abs(x - (xOffset + newWidth / 2));
					double distY = Math.abs(y - (yOffset + newHeight / 2));
					double distance = Math.sqrt(distX * distX + distY * distY);
					
					// If the distance is greater than the stretch radius, we stretch the edge pixels
					if (distance > stretchRadius) {
						// Get the closest edge pixel coordinates
						int srcX = Math.max(0, Math.min(newWidth - 1, x - xOffset));
						int srcY = Math.max(0, Math.min(newHeight - 1, y - yOffset));
						
						// Copy the edge pixel from the scaled image to the letterboxed area
						imgOut.setRGB(x, y, scaledImage.getRGB(srcX, srcY));
					}
				}
			}
		}
		
		g2.dispose();
	}

	private void stretchEdges(Graphics2D g2, BufferedImage imgIn, int newWidth, int newHeight, int xOffset, int yOffset, int outWidth, int outHeight, double stretchRadius) {
		int centerX = xOffset + newWidth / 2;
		int centerY = yOffset + newHeight / 2;
	
		// Stretch horizontally
		for (int x = 0; x < outWidth; x++) {
			if (x < xOffset || x >= xOffset + newWidth) {
				// Pixels outside the original image (left and right letterbox)
				double distanceX = (x < xOffset) ? xOffset - x : x - (xOffset + newWidth);
				double stretchX = distanceX / stretchRadius;
	
				for (int y = 0; y < newHeight; y++) {
					int srcX = (x < xOffset) ? 0 : newWidth - 1;
					int srcY = y;
					int targetX = (x < xOffset) ? x : outWidth - (xOffset + newWidth) + xOffset + x;
					g2.drawImage(imgIn.getSubimage(srcX, srcY, 1, 1), targetX, yOffset + y, 1, 1, null);
				}
			}
		}
		
		// Stretch vertically
		for (int y = 0; y < outHeight; y++) {
			if (y < yOffset || y >= yOffset + newHeight) {
				// Pixels outside the original image (top and bottom letterbox)
				double distanceY = (y < yOffset) ? yOffset - y : y - (yOffset + newHeight);
				double stretchY = distanceY / stretchRadius;
	
				for (int x = 0; x < newWidth; x++) {
					int srcX = x;
					int srcY = (y < yOffset) ? 0 : newHeight - 1;
					int targetY = (y < yOffset) ? y : outHeight - (yOffset + newHeight) + yOffset + y;
					g2.drawImage(imgIn.getSubimage(srcX, srcY, 1, 1), xOffset + x, targetY, 1, 1, null);
				}
			}
		}
	}
	

//x = ((int)nonLinearStep(width, outWidth, x+1));
	private void downSampleNonLinear(int width, int height,int outWidth, int outHeight, BufferedImage imgIn, BufferedImage imgOut)
	{
		//int y = 0;
		int outY = 0;
		while(outY < outHeight)
		{
			//int x = 0;
			int outX = 0;
			while(outX < outWidth)
			{
				//System.out.println("herex " +x);
				// int color = imgIn.getRGB(clamp((nonLinearStepFixed(width/2, outWidth/2, outX-(outWidth/2)))+(width/2), 0, width-1),clamp((nonLinearStepFixed(height/2, outHeight/2, outX-(outWidth/2)))+(height/2),0,height-1));
				// int color = imgIn.getRGB(clamp((nonLinearStepFixed(width/2, outWidth/2, outX-(outWidth/2)))+(width/2), 0, width-1),clamp((nonLinearStepFixed(height/2, outHeight/2, outY-(outHeight/2)))+(height/2),0,height-1));

				// Apply the reversed non-linear transformation

				//Works for pos x and negative y
				//int newX = clamp(linearStep(width/2, outWidth/2, (outX-(outWidth/2)))+(width/2), 0, width - 1);
				int newX = 0;
				if((outX-(outWidth/2)) >= 0){
					newX = clamp(nonLinearStepFixed2(width/2, height/2, outWidth/2, (outX-(outWidth/2)))+(width/2), 0, width - 1);
				} else {
					newX = clamp(-nonLinearStepFixed2(width/2, height/2, outWidth/2, -(outX-(outWidth/2)))+(width/2), 0, width - 1);
				}

				int newY = 0;
				if((outY-(outHeight/2)) >= 0){
					newY = clamp(nonLinearStepFixed2(height/2, height/2, outHeight/2, (outY-(outHeight/2)))+(height/2), 0, height - 1);
				} else {
					newY = clamp(-nonLinearStepFixed2(height/2,height/2, outHeight/2, -(outY-(outHeight/2)))+(height/2), 0, height - 1);;
				}
				//int newX = clamp(nonLinearStepFixed(width/2, outWidth/2, (outX-(outWidth/2)))+(width/2), 0, width - 1);
				//int newY = clamp(linearStep(height/2, outHeight/2, (outY-(outHeight/2)))+(height/2), 0, height - 1);
	
				// Get the color from the original image
				int color = imgIn.getRGB(newX, newY);
				// Components will be in the range of 0..255:
				byte r = (byte)((color & 0xff0000) >> 16);
				byte g = (byte)((color & 0xff00) >> 8);
				byte b = (byte)(color & 0xff);
				
				byte a = (byte)((color & 0xff000000) >> 24);

				//int pixOut = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				int pix = ((a << 24) + (r << 16) + (g << 8) + b);
				//int pix = 0xff32A852;
				imgOut.setRGB(outX,outY,pix);
				// x = ((int)nonLinearStep(width, outWidth, x+1));
				//x = (linearStep(width, outWidth, x+1));
				outX++;
			}
			//System.out.println("herey " +y);
			// y = ((int)nonLinearStep(height, outHeight, y+1));
			//y = (linearStep(height, outHeight, y+1));
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
		int width = 4000;
		int height = 3000;

		// Read in the specified image
		BufferedImage res = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		readImageRGB(width, height, args[0], res);

		// //Downsample1
		// int outWidth = 640;
		// int outHeight = 400;
		// BufferedImage ds1Out = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		// downSample1(width, height, outWidth, outHeight, res, ds1Out);
		// res = ds1Out;

		//Downsample2
		//GaussBlur
		// double rho = 3;

		// BufferedImage imgBlurOut = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		// gaussBlur(width, height, (int)(3*rho), rho, res, imgBlurOut);
		// res = imgBlurOut;
		
		// // //Then Ds1
		// int outWidth = 640;
		// int outHeight = 400;
		// BufferedImage ds1Out = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		// downSample1(width, height, outWidth, outHeight, res, ds1Out);
		// res = ds1Out;

		// //Upsample
		// int outWidth = 1920;
		// int outHeight = 1080;
		// BufferedImage us1Out = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		// upSample1(width, height, outWidth, outHeight, res, us1Out);
		// res = us1Out;

		// //Upsample
		// int outWidth = 1920;
		// int outHeight = 1080;
		// BufferedImage biOut = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		// bilinear(width, height, outWidth, outHeight, res, biOut);
		// res = biOut;

		//PAR
		int outWidth = 1920;
		int outHeight = 1080;
		BufferedImage ds1Out = new BufferedImage(outWidth, outHeight, BufferedImage.TYPE_INT_RGB);
		//downSampleFUCK(width, height, outWidth, outHeight, res, ds1Out, 400);
		downSampleFirst(width, height, outWidth, outHeight, res, ds1Out);
		// Graphics2D sGraphics2d = new G
		// stretchEdges(null, ds1Out, outWidth, outHeight, width, height, outWidth, outHeight, outHeight);
		res = ds1Out;

		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		// lbIm1 = new JLabel(new ImageIcon(imgOne));
		// lbIm1 = new JLabel(new ImageIcon(imgDS1Out));
		lbIm1 = new JLabel(new ImageIcon(res));

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
			ImageIO.write(res, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// try {
		// 	ImageIO.write(imgBlur, "png", outputfile);
		// } catch (IOException e) {
		// 	// TODO Auto-generated catch block
		// 	e.printStackTrace();
		// }
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

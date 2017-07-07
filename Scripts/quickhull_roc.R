
## 2016, Andrej Aderhold
##
## Quickhull algorithm adapted for ROC curve
##
## This algorithm finds the convex hull of the ROC curve, for points
## located at the upper left triangle of the curve.
##
## Arguments:
##        x   : FPR (False Positive Rate) = 1 - Specificity
##        y   : Sensitivity
##


quickhull_roc <- function(x,y) {

    ## convert FPR into Specificity
    x = 1 - x
    
    ## set to TRUE, will plot the hull segments and points selected for the hull
    DEBUG_OUTPUT = FALSE
    
    ## initial coordinates that divide the ROC area
    hull_x_set = c(1,0)
    hull_y_set = c(0,1)

    ## flag to remain in the main while loop as long as new hull coordinates are found
    hull_points_added = TRUE


    ## main while loop
    ##
    ## In each iteration of this loop, all segments of the hull are checked for distant points
    ## following the divide and conquer algorithm of Quickhull.
    ## If new points are found, the hull is updated and the whole search is repeated.
    ##
    ## Note: to speed up things I could flag hull segments that have been already checked as done.
    ##
    while(hull_points_added) {

        ## Reset this flag. Whenever a new hull point was found it is set to TRUE below.
        ## If none was found, the while loop finishes.
        hull_points_added = FALSE


        ## save to buffer, new points go into newset and are written back at the end of this loop
        hull_x_newset = hull_x_set
        hull_y_newset = hull_y_set
   

        ## save the points that were added to the hull into here, these are deleted right after this loop
        points_to_del = c()
    
        ## Try to find new hull points along the set of points in hull_x_set and hull_y_set
        for(mj in 1:(length(hull_x_set) - 1)) {

            ## extract the pair of adjacent hull coordinates, this is one segment of the hull
            hull_x = c(hull_x_set[mj], hull_x_set[mj+1])
            hull_y = c(hull_y_set[mj], hull_y_set[mj+1])

            ## We need to find out if a point of interest is outside the hull or not.
            ## If it is outside, it could be a potential hull candidate. 
            ## This filter step and involves the calculation of the hull segment line, that connects the
            ## current two points of the hull.
            ##
            ## The slope : delta_y / delta_x
            m = (hull_y[2] - hull_y[1]) / (hull_x[2] - hull_x[1])

            ## the intercept given one of the hull points, i.e.  y - (m*x) = b , where b is the intercept
            intercept = hull_y[1] - ( m * hull_x[1] )

            ## remember the most distant points from the hull segment 
            distmax = 0
            max_i = -1

            ## loop over each potential hull point
            ##  1. check if it outside the hull, if not, skip this point
            ##  2. if outside, calculate the distance to both hull points
            ##  3. remember the point with the furthest distance and add to hull
            ##
            for ( i in 1:length(x)) {


                ## find out if point is equal to boundary points of hull segment
                ## skip this if applies
                if((x[i] == hull_x[1]) && (y[i] == hull_y[1])) ## check left boundary
                    next
                if((x[i] == hull_x[2]) && (y[i] == hull_y[2])) ## check right boundary
                    next
                
                
                ## find out if point (x[i],y[i]) is above the hull segment
                y_star = intercept + m * x[i]
                
                ## if the line is vertical, detect it and look if point has larger x value 
                if (hull_x[1] == hull_x[2]) {  ## slope is Inf 
                                       
                    ## if point in question is smaller than the current hull fragment, skip it
                    if (hull_x[1] >= x[i])
                        next
                
                } else if (y_star >= y[i])  ## if point has smaller y value, skip to next point
                    next
            
                ## if point was outside hull, calculate the euclidian distance to both hull points
                dist1 = sqrt( ( x[i] - hull_x[1])^2 + (y[i] - hull_y[1])^2 )
                dist2 = sqrt( ( x[i] - hull_x[2])^2 + (y[i] - hull_y[2])^2 )

                ## sum distances
                disttotal = dist1 + dist2

                ##
                ## Do some plotting of points and hull segment to see whats going on
                ##
                if (DEBUG_OUTPUT) {

                    ## plot all points
                    plot(x, y, xlim=c(1,0), ylim=c(0,1))
                    title(sprintf("distance for point %i: %g", i, disttotal))
                    #par(new=T)
                    
                    ## plot point with maximum distance if found
                    if(max_i > 0) {
                        par(new=T)
                        plot(x[max_i], y[max_i], pch=2, lwd=2, col="blue", xlim=c(1,0), ylim=c(0,1))
                        
                    }
                    
                    ## plot current point
                    par(new=T)
                    plot(x[i]    , y[i]    , pch=1, lwd=2, col="green", xlim=c(1,0), ylim=c(0,1))
                    
                    ## plot the hull segment, i.e. the divide line, used to determine distant point
                    par(new=T)
                    plot(hull_x, hull_y, type="b", xlim=c(1,0), ylim=c(0,1))

                    ## plot the full hull in grey as reference
                    par(new=T)
                    plot(hull_y_newset ~ hull_x_newset, col = "grey",type = "b", lwd = 1, lty = 2, xlim = c(1,0), ylim = c(0,1))

                }


                
                ## .. and check if it is the maximum distance
                if (distmax < disttotal) {
                    ## if yes, save the maximum distance and the index of this point
                    distmax = disttotal
                    max_i = i
                }

                                    
            } ## end loop check points if element of hull


            ##
            ## if a distant point was found, insert it into the hull
            ##
            if (max_i > 0) {
                
                hull_x_newset = append(hull_x_newset, x[max_i], mj)
                hull_y_newset = append(hull_y_newset, y[max_i], mj)

                if (DEBUG_OUTPUT) {
                    par(new=T)
                    plot(x[max_i], y[max_i], pch=3, lwd=5, col="red", xlim=c(1,0), ylim=c(0,1))

                    ## plot the full hull in grey as reference
                    par(new=T)
                    plot(hull_y_newset ~ hull_x_newset, col = "grey",type = "b", lwd = 1, lty = 2, xlim = c(1,0), ylim = c(0,1))

                }
                    
                
                
                ## remove this hull point from point set
                ## only remember the indices and do the deletion later, otherwise we
                ## mess up the loop of the hull points in which we are still in. 
                points_to_del = c(points_to_del, max_i)
                
                ## flag that we found a new point so the while loop executes again
                hull_points_added = TRUE

                ## Break, so that we can replace the last hull with the new hull
                ## If this does not happen, we get troubles with inserting new points
                ## into hull_*_newset using the index 'mj' from hull_*_set.
                break
            } 

            if(DEBUG_OUTPUT) {
                print.table(hull_x_newset)
                print.table(hull_y_newset)
            }
                            
        } ## end loop over hull points

        ## now we can update the hull with the points we previously added to 'hull_[x|y]_set'
        hull_x_set = hull_x_newset
        hull_y_set = hull_y_newset

        ## take out points that were previously added to the hull 
        if (length(points_to_del) > 0) {
            x = x[- points_to_del]
            y = y[- points_to_del]
        }
        
    } ## end of main while loop

    cat(" * Convex hull calculated.\n")

    ## convert the set of x values back from specificity to FPR
    fpr_hull = 1 - hull_x_set
    
    return(list(x = fpr_hull, y = hull_y_set))
}

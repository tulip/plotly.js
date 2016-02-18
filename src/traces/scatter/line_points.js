/**
* Copyright 2012-2016, Plotly, Inc.
* All rights reserved.
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


'use strict';

var Axes = require('../../plots/cartesian/axes');


module.exports = function linePoints(d, opts) {
    var xa = opts.xaxis,
        ya = opts.yaxis,
        connectGaps = opts.connectGaps,
        baseTolerance = opts.baseTolerance,
        linear = opts.linear,
        segments = [],
        badnum = Axes.BADNUM,
        minTolerance = 0.2, // fraction of tolerance "so close we don't even consider it a new point"
        pts = new Array(d.length),
        pti = 0,
        i,

        // pt variables are pixel coordinates [x,y] of one point
        clusterStartPt, // these four are the outputs of clustering on a line
        clusterEndPt,
        clusterHighPt,
        clusterLowPt,
        thisPt, // "this" is the next point we're considering adding to the cluster

        clusterRefDist,
        clusterHighFirst, // did we encounter the high point first, then a low point, or vice versa?
        clusterUnitVector, // the first two points in the cluster determine its unit vector
                           // so the second is always in the "High" direction
        thisVector, // the pixel delta from clusterStartPt

        // val variables are (signed) pixel distances along the cluster vector
        clusterHighVal,
        clusterLowVal,
        thisVal,

        // deviation variables are (signed) pixel distances normal to the cluster vector
        clusterMinDeviation,
        clusterMaxDeviation,
        thisDeviation;

    // turn one calcdata point into pixel coordinates
    function getPt(index) {
        var x = xa.c2p(d[index].x),
            y = ya.c2p(d[index].y);
        if(x === badnum || y === badnum) return false;
        return [x, y];
    }

    // if we're off-screen, increase tolerance over baseTolerance
    function getTolerance(pt) {
        var xFrac = pt[0] / xa._length,
            yFrac = pt[1] / ya._length;
        return (1 + 10 * Math.max(0, -xFrac, xFrac - 1, -yFrac, yFrac - 1)) * baseTolerance;
    }

    function ptDist(pt1, pt2) {
        var dx = pt1[0] - pt2[0],
            dy = pt1[1] - pt2[1];
        return Math.sqrt(dx * dx + dy * dy);
    }

    // loop over ALL points in this trace
    for(i = 0; i < d.length; i++) {
        clusterStartPt = getPt(i);
        if(!clusterStartPt) continue;

        pti = 0;
        pts[pti++] = clusterStartPt;

        // loop over one segment of the trace
        for(i++; i < d.length; i++) {
            clusterHighPt = getPt(i);
            if(!clusterHighPt) {
                if(connectGaps) continue;
                else break;
            }

            // can't decimate if nonlinear line shape
            // TODO: we *could* decimate [hv]{2,3} shapes if we restricted clusters to horz or vert again
            // but spline would be verrry awkward to decimate
            if(!linear) {
                pts[pti++] = clusterHighPt;
                continue;
            }

            clusterRefDist = ptDist(clusterHighPt, clusterStartPt);

            if(clusterRefDist < getTolerance(clusterHighPt) * minTolerance) continue;

            clusterUnitVector = [
                (clusterHighPt[0] - clusterStartPt[0]) / clusterRefDist,
                (clusterHighPt[1] - clusterStartPt[1]) / clusterRefDist
            ];

            clusterLowPt = clusterStartPt;
            clusterHighVal = clusterRefDist;
            clusterLowVal = clusterMinDeviation = clusterMaxDeviation = 0;
            clusterHighFirst = false;
            clusterEndPt = clusterHighPt;

            // loop over one cluster of points that collapse onto one line
            for(i++; i < d.length; i++) {
                thisPt = getPt(i);
                if(!thisPt) {
                    if(connectGaps) continue;
                    else break;
                }
                thisVector = [
                    thisPt[0] - clusterStartPt[0],
                    thisPt[1] - clusterStartPt[1]
                ];
                // cross product (or dot with normal to the cluster vector)
                thisDeviation = thisVector[0] * clusterUnitVector[1] - thisVector[1] * clusterUnitVector[0];
                clusterMinDeviation = Math.min(clusterMinDeviation, thisDeviation);
                clusterMaxDeviation = Math.max(clusterMaxDeviation, thisDeviation);

                if(clusterMaxDeviation - clusterMinDeviation > getTolerance(thisPt)) break;

                clusterEndPt = thisPt;
                thisVal = thisVector[0] * clusterUnitVector[0] + thisVector[1] * clusterUnitVector[1];

                if(thisVal > clusterHighVal) {
                    clusterHighVal = thisVal;
                    clusterHighPt = thisPt;
                    clusterHighFirst = false;
                } else if(thisVal < clusterLowVal) {
                    clusterLowVal = thisVal;
                    clusterLowPt = thisPt;
                    clusterHighFirst = true;
                }
            }

            // insert this cluster into pts
            // we've already inserted the start pt, now check if we have high and low pts
            if(clusterHighFirst) {
                pts[pti++] = clusterHighPt;
                if(clusterEndPt !== clusterLowPt) pts[pti++] = clusterLowPt;
            } else {
                if(clusterLowPt !== clusterStartPt) pts[pti++] = clusterLowPt;
                if(clusterEndPt !== clusterHighPt) pts[pti++] = clusterHighPt;
            }
            // and finally insert the end pt
            pts[pti++] = clusterEndPt;

            // have we reached the end of this segment?
            if(i >= d.length || !thisPt) break;

            // otherwise we have an out-of-cluster point to insert as next clusterStartPt
            pts[pti++] = thisPt;
            clusterStartPt = thisPt;
        }

        segments.push(pts.slice(0, pti));
    }

    console.log('segments', segments);
    // console.log('xaxis', opts.xaxis.range);
    // console.log('yaxis', opts.yaxis.range);
    test(d, opts);

    return segments;
};


function test(points, opts) {
    var xRange = opts.xaxis.range,
        yRange = opts.yaxis.range;

    // Plot bounding points.
    var bounds = [
        { x: xRange[0], y: yRange[0] },
        { x: xRange[0], y: yRange[1] },
        { x: xRange[1], y: yRange[1] },
        { x: xRange[1], y: yRange[0] }
    ];

    var segmentss = [],
        added = false;

    for(var i = 0; i < points.length; i++){
        var point = points[i];

        if(insideBox(point) && !added) segmentss.push(point);
        if(i + 1 >= points.length) break;

        var line = [points[i], points[i + 1]],
            lineSlope = (line[1].y - line[0].y) / (line[1].x - line[0].x),
            slopes = [];

        for(var k = 0; k < bounds.length; k++){
            var slope = (bounds[k].y - line[0].y) / (bounds[k].x - line[0].x);
            slopes.push(slope);
        }

        slopes.sort(function(a, b){
            return a - b;
        });

        var sweepsPlotArea = lineSlope >= slopes[0] && lineSlope <= slopes[slopes.length - 1];

        if(sweepsPlotArea || insideBox(line[1])){
            if(!added){
                segmentss.push(points[i]);
            }
            segmentss.push(points[i + 1]);
            added = true;
        }else{
            added = false;
        }
    }

    function insideBox(p){
        return (
            bounds[0].x < p.x && p.x < bounds[2].x &&
            bounds[0].y < p.y && p.y < bounds[2].y
        );
    }

    function getPt(p) {
        var x = opts.xaxis.c2p(p.x),
            y = opts.yaxis.c2p(p.y);
        if(x === Axes.badnum || y === Axes.badnum) return false;
        return [x, y];
    }

    // Decimate now.
    for(var j = 0; j < segmentss.length; j++){
        segmentss[j] = getPt(segmentss[j]);
    }

    console.log('my segments', [segmentss]);

    return [segmentss.reverse()];
}

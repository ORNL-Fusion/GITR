FN0 = 50;
side_base=0.03;
angle = 30;
sub_width = 2*side_base;
slope = tan(angle);
y = 1/sqrt(slope*slope+1);
x = -slope/sqrt(slope*slope+1);
difference (){
    translate([0, 0, 0])
    cube([side_base,side_base,2*side_base], center = true,$fn=FN0);

translate([-0.5*x*sub_width, 0,-0.5*y*sub_width]){
rotate(a=90-angle,v=[0,1,0]) { 
    
    cube([sub_width,sub_width,2*sub_width], center = true,$fn=FN0);
    }

}

}

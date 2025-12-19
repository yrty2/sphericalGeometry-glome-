let moveVector=[0,0,0];
let rotation=[0,0,0];
const r=300;
let frames=0;
let auto=false;
const camera=[0,0,0];//pos-cam->viewpos
const canvas=document.querySelector(".canvas");
canvas.width=screen.width-4;
canvas.height=screen.height-4;
canvas.style.border="2px solid";
const ctx=canvas.getContext("2d");
function clip(u,b){
    const v=projection.stereographic3D(u);
    v.x-=camera[0];
    v.y=-(v.y+camera[1]);
    v.z-=camera[2];
    if(v.z>0){
        const p=projection.perspective(v).scale(canvas.height).add(new cartesian2D(canvas.width/2,canvas.height/2));
        if(b){
            return new cartesian(p.x,p.y,canvas.height/v.z);
        }
    return p;
    }
}
function point(p){
    if(p!==undefined){
    ctx.beginPath();
    ctx.arc(p.x,p.y,5,0,2*Math.PI);
    ctx.fill();
    ctx.closePath();
    }
}
function pointwithSize(p){
    if(p!==undefined && p.z<200){
    ctx.beginPath();
    ctx.arc(p.x,p.y,p.z*4,0,2*Math.PI);
    ctx.fill();
    ctx.closePath();
    }
}
function plot(v){
    pointwithSize(clip(v,true));
}
const points=[];
const colors=[];
function translate(){
    ctx.fillStyle="#000000";
    ctx.fillRect(0,0,canvas.width,canvas.height);
    for(let k=0; k<points.length; ++k){
        let p=points[k];
        p.rotate(rotation[0],rotation[1],rotation[2]);
        p.translate(moveVector[0],moveVector[1],moveVector[2]);
        ctx.fillStyle=`rgba(${255*colors[k][0]},${255*colors[k][1]},${255*colors[k][2]},0.66)`;
        plot(p);
        //cube(p,10);
    }
    keycontrol();
    frames++;
    requestAnimationFrame(translate);
}
translate();
generate();
function generate(){
    const S=30;//12
    for(let x=0; x<S; ++x){
    for(let y=0; y<S; ++y){
    for(let z=0; z<S; ++z){
        colors.push([x/S,y/S,z/S]);
        points.push(new spherical4D(r,new cartesian4D(
            Math.cos(2*Math.PI*x/S)*Math.sin(2*Math.PI*y/S)*Math.sin(2*Math.PI*z/S),
            Math.sin(2*Math.PI*x/S)*Math.sin(2*Math.PI*y/S)*Math.sin(2*Math.PI*z/S),
            Math.cos(2*Math.PI*y/S)*Math.sin(2*Math.PI*z/S),
            Math.cos(2*Math.PI*z/S)).scale(r)));
    }
    }
    }
}
function line(a,dir,length){
    //dirは３次元ベクトル
    //aからbへの直線
    const d=100;
    const now=a.copy();
    const s=d*Math.sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
    ctx.beginPath();
    for(let k=0; k<d; ++k){
        const p=clip(now);
        if(p){
        ctx.lineTo(p.x,p.y);
        }
        now.translate(length*dir[0]/s,length*dir[1]/s,length*dir[2]/s);
    }
    ctx.stroke();
    ctx.closePath();
}
function cube(p,s){
    line(p.translateBack(-s/2,-s/2,-s/2),[1,0,0],s);
    line(p.translateBack(s/2,-s/2,-s/2),[0,1,0],s);
    line(p.translateBack(s/2,s/2,-s/2),[-1,0,0],s);
    line(p.translateBack(-s/2,s/2,-s/2),[0,-1,0],s);

    line(p.translateBack(-s/2,-s/2,s/2),[1,0,0],s);
    line(p.translateBack(s/2,-s/2,s/2),[0,1,0],s);
    line(p.translateBack(s/2,s/2,s/2),[-1,0,0],s);
    line(p.translateBack(-s/2,s/2,s/2),[0,-1,0],s);

    line(p.translateBack(-s/2,-s/2,-s/2),[0,0,1],s);
    line(p.translateBack(s/2,-s/2,-s/2),[0,1,1],s);
    line(p.translateBack(s/2,s/2,-s/2),[0,0,1],s);
    line(p.translateBack(-s/2,s/2,-s/2),[0,0,1],s);
}

"use strict";
function rgb(h,s,l){h=((h%360)+360)%360;s/=100;l/=100;const a=s*Math.min(l,1-l);const f=n=>{const k=(n+h/30)%12;return Math.round(255*(l-a*Math.max(-1,Math.min(k-3,9-k,1))))};return[f(0),f(8),f(4)]}
const themes={aurora:[194,30,152,279,12,74,327,230],ember:[20,46,350,280,64,8,200,322],glacier:[180,203,233,269,156,330,52,119],spectrum:[10,65,135,185,225,290,338,102]};
self.onmessage=({data})=>{const {id,co,roots,view,width,height,palette,maxIter}=data;try{const hues=themes[palette]||themes.aurora;const colors=roots.map((_,j)=>Array.from({length:maxIter+1},(_,k)=>rgb(hues[j%hues.length]+Math.floor(j/hues.length)*17,77,Math.max(14,75-54*Math.pow(k/maxIter,.67)))));const px=new Uint8ClampedArray(width*height*4);const R=view.R,asp=height/width;
for(let row=0;row<height;row++){const im=view.y+(1-2*(row+.5)/height)*R*asp;for(let col=0;col<width;col++){let zr=view.x+(2*(col+.5)/width-1)*R,zi=im,match=-1,k=0;
for(;k<maxIter;k++){let pr=co[0],pi=0,dr=0,di=0;for(let a=1;a<co.length;a++){const ndr=dr*zr-di*zi+pr,ndi=dr*zi+di*zr+pi;const npr=pr*zr-pi*zi+co[a],npi=pr*zi+pi*zr;dr=ndr;di=ndi;pr=npr;pi=npi}
const denom=dr*dr+di*di;if(!Number.isFinite(denom)||denom<1e-28)break;const sr=(pr*dr+pi*di)/denom,si=(pi*dr-pr*di)/denom;zr-=sr;zi-=si;if(!Number.isFinite(zr)||!Number.isFinite(zi)||Math.abs(zr)+Math.abs(zi)>1e13)break;
if(Math.hypot(sr,si)<1e-8||k>2){for(let j=0;j<roots.length;j++)if(Math.hypot(zr-roots[j][0],zi-roots[j][1])<2e-5){match=j;break}if(match>=0)break}
}const i=(row*width+col)*4;if(match>=0){const a=colors[match][Math.min(maxIter,k)];px[i]=a[0];px[i+1]=a[1];px[i+2]=a[2]}else{px[i]=7;px[i+1]=13;px[i+2]=24}px[i+3]=255}}
self.postMessage({id,width,height,pixels:px.buffer},[px.buffer])}catch(e){self.postMessage({id,error:String(e.message||e)})}};
"use strict";
const DATES=[
["25/12/1642","Isaac Newton · Julian","Mathematician"],["04/01/1643","Isaac Newton · Gregorian","Mathematician"],["15/04/1707","Leonhard Euler","Mathematician"],["30/04/1777","Carl Friedrich Gauss","Mathematician"],["16/08/1821","Arthur Cayley","Mathematician"],["29/04/1854","Henri Poincaré","Mathematician"],["28/02/1878","Pierre Fatou","Mathematician"],["03/02/1893","Gaston Julia","Mathematician"],["20/11/1924","Benoit Mandelbrot","Mathematician"],["05/07/1687","Principia · Julian","Scientific milestone"],["15/07/1687","Principia · Gregorian","Scientific milestone"],["07/09/1822","Independence of Brazil","Historical milestone"],["23/09/1846","Observation of Neptune","Scientific milestone"],["13/05/1888","Abolition of slavery in Brazil","Historical milestone"],["21/01/1889","King Oscar II's prize","Mathematical history"],["15/11/1889","Proclamation of Brazil's Republic","Historical milestone"],["30/11/1889","Poincaré stops printing","Mathematical history"],["08/08/1900","Hilbert's problems","Mathematical history"],["23/10/1906","Flight of the 14-bis","Scientific milestone"],["29/05/1919","Sobral eclipse","Scientific milestone"],["15/10/1952","Founding of IMPA","Mathematical history"],["29/12/1972","Lorenz's butterfly talk","Scientific milestone"],["16/05/1980","Founding of LNCC","Mathematical history"],["05/10/1988","Brazilian Constitution","Historical milestone"],["23/06/1993","Wiles's announcement","Mathematical history"],["22/04/1500","Arrival of Cabral · Julian","Historical milestone"],["02/05/1500","Arrival of Cabral · Gregorian","Historical milestone"],["04/10/1582","Last Julian calendar day","Historical milestone"],["15/10/1582","First Gregorian calendar day","Historical milestone"],["20/07/1969","Apollo 11 · Brasília","Scientific milestone"],["21/07/1969","Apollo 11 · UTC","Scientific milestone"],["14/09/2015","Gravitational waves · detection","Scientific milestone"],["11/02/2016","Gravitational waves · announcement","Scientific milestone"],["01/01/2000","New Year's Day 2000","Historical milestone"],["30/08/1984","Personal date · author's birth","Personal"],["03/07/1986","Personal date","Personal"],["27/04/1952","Personal date","Personal"],["15/08/1950","Personal date","Personal"]];
const $=id=>document.getElementById(id);let worker=new Worker("worker.js");const callbacks=new Map();let nextId=1;let current=null;let visibleDates=12;
worker.onmessage=({data})=>{const fn=callbacks.get(data.id);if(fn){callbacks.delete(data.id);fn(data)}};
worker.onerror=()=>{$("renderStatus").textContent="Renderer error";};
function isValidDate(s){let m=/^(\d{2})\/(\d{2})\/(\d{4})$/.exec(s.trim());if(!m)return false;let d=+m[1],mo=+m[2],y=+m[3];if(y<1||mo<1||mo>12)return false;let leap=y%4===0&&(y%100!==0||y%400===0);return d>=1&&d<=[31,leap?29:28,31,30,31,30,31,31,30,31,30,31][mo-1]}
function coeffs(s){if(!isValidDate(s))throw Error("Enter a valid date as DD/MM/YYYY.");let a=[...s.replaceAll("/","")].map(Number);while(a[0]===0)a.shift();return a}
function polyHTML(co){const deg=co.length-1;const terms=[];co.forEach((a,i)=>{if(!a)return;const k=deg-i;terms.push((a===1&&k>0?"":String(a))+(k===0?"":k===1?"z":"z<sup>"+k+"</sup>"))});return "p(z) = "+terms.join(" + ")}
function add(a,b){return[a[0]+b[0],a[1]+b[1]]}function mul(a,b){return[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]]}
function div(a,b){let q=b[0]*b[0]+b[1]*b[1];return [(a[0]*b[0]+a[1]*b[1])/q,(a[1]*b[0]-a[0]*b[1])/q]}
function rootsOf(co){let q=co.slice(),zero=0;while(q.length>2&&q.at(-1)===0){q.pop();zero++}const n=q.length-1;const c=q.map(v=>v/q[0]);const rad=Math.max(1.25,Math.min(12,1+Math.max(...c.slice(1).map(Math.abs))));let z=Array.from({length:n},(_,i)=>[rad*Math.cos(2*Math.PI*(i+.137)/n),rad*Math.sin(2*Math.PI*(i+.137)/n)]);
for(let k=0;k<400;k++){let zs=[],deltaMax=0;for(let i=0;i<n;i++){let p=[c[0],0],product=[1,0];for(let j=1;j<c.length;j++)p=add(mul(p,z[i]),[c[j],0]);for(let j=0;j<n;j++)if(j!==i)product=mul(product,[z[i][0]-z[j][0],z[i][1]-z[j][1]]);if(Math.hypot(...product)<1e-22)throw Error("Ill-conditioned root calculation for this date.");const step=div(p,product);zs.push([z[i][0]-step[0],z[i][1]-step[1]]);deltaMax=Math.max(deltaMax,Math.hypot(...step))}z=zs;if(deltaMax<1e-10)break}
if(z.some(v=>!Number.isFinite(v[0])||!Number.isFinite(v[1])))throw Error("Root approximation failed.");
const unique=[];for(const root of z)if(!unique.some(r=>Math.hypot(r[0]-root[0],r[1]-root[1])<1e-5))unique.push(root);if(zero)unique.push([0,0]);return unique}
function bounds(roots){const re=roots.map(r=>r[0]).sort((a,b)=>a-b),im=roots.map(r=>r[1]).sort((a,b)=>a-b);const cx=re[Math.floor(re.length/2)],cy=im[Math.floor(im.length/2)];let d=roots.map(r=>Math.hypot(r[0]-cx,r[1]-cy)).sort((a,b)=>a-b);const med=d[Math.floor(d.length/2)];let kept=roots.filter(r=>Math.hypot(r[0]-cx,r[1]-cy)<3*Math.max(.25,med));if(kept.length<2)kept=roots;let x=kept.reduce((v,r)=>v+r[0],0)/kept.length,y=kept.reduce((v,r)=>v+r[1],0)/kept.length,R=Math.max(...kept.map(r=>Math.hypot(r[0]-x,r[1]-y)));return {x,y,R:Math.max(.35,1.4*R)}}
function render(opts){return new Promise(resolve=>{const id=nextId++;callbacks.set(id,resolve);worker.postMessage({id,...opts})})}
function paint(canvas,data,roots,view){let ctx=canvas.getContext("2d");canvas.width=data.width;canvas.height=data.height;ctx.putImageData(new ImageData(new Uint8ClampedArray(data.pixels),data.width,data.height),0,0);ctx.strokeStyle="#fff";ctx.fillStyle="#fff";ctx.lineWidth=Math.max(1,canvas.width/550);roots.forEach(([x,y])=>{const px=(x-view.x)/view.R*canvas.width/2+canvas.width/2,py=(view.y-y)/view.R*canvas.width/2+canvas.height/2;if(px>=0&&px<canvas.width&&py>=0&&py<canvas.height){ctx.beginPath();ctx.arc(px,py,Math.max(1.4,canvas.width/280),0,Math.PI*2);ctx.stroke()}})}
const params=new URLSearchParams(location.search);if(isValidDate(params.get("date")||""))$("date").value=params.get("date");if(["aurora","ember","glacier","spectrum"].includes(params.get("palette")))$("palette").value=params.get("palette");
let run=0;
async function generate(newDate=true){const serial=++run;const date=$("date").value.trim();$("dateError").textContent="";let co,roots;try{co=coeffs(date);roots=rootsOf(co)}catch(e){$("dateError").textContent=e.message;return}let view=newDate||!current?bounds(roots):current.view;current={date,co,roots,view};$("polynomial").innerHTML=polyHTML(co);$("degree").textContent=co.length-1;$("rootCount").textContent=roots.length;$("windowLabel").textContent="±"+view.R.toFixed(4);$("renderStatus").textContent="Computing…";
const size=+$("resolution").value;const result=await render({co,roots,view,width:size,height:size,palette:$("palette").value,maxIter:85});if(serial!==run)return;if(result.error){$("renderStatus").textContent=result.error;return}paint($("fractalCanvas"),result,roots,view);$("renderStatus").textContent="Image ready · "+size+" × "+size;const u=new URL(location.href);u.searchParams.set("date",date);u.searchParams.set("palette",$("palette").value);history.replaceState(null,"",u)}
$("fractalForm").onsubmit=e=>{e.preventDefault();generate(true)};$("heroForm").onsubmit=e=>{e.preventDefault();$("date").value=$("heroDate").value;document.querySelector("#studio").scrollIntoView({behavior:"smooth"});generate(true)};
document.querySelectorAll("[data-preset]").forEach(b=>b.onclick=()=>{$("date").value=b.dataset.preset;generate(true)});
function zoom(factor,pt=null){if(!current)return;const v=current.view;if(pt){v.x+=(pt[0]*2-1)*v.R;v.y-=(pt[1]*2-1)*v.R}v.R*=factor;generate(false)}
$("zoomIn").onclick=()=>zoom(.5);$("zoomOut").onclick=()=>zoom(2);$("resetView").onclick=()=>{if(current){current.view=bounds(current.roots);generate(false)}};
$("fractalCanvas").addEventListener("click",e=>{let b=e.currentTarget.getBoundingClientRect();zoom(.5,[(e.clientX-b.left)/b.width,(e.clientY-b.top)/b.height])});
$("fractalCanvas").addEventListener("keydown",e=>{if(e.key==="Enter"||e.key==="+"){zoom(.5);e.preventDefault()}});
$("download").onclick=()=>{if(!current)return;const a=document.createElement("a");a.download="BirthdayFractal-"+current.date.replaceAll("/","-")+".png";a.href=$("fractalCanvas").toDataURL("image/png");a.click()};
$("copyLink").onclick=async()=>{const b=$("copyLink");try{await navigator.clipboard.writeText(location.href);b.textContent="Copied!";setTimeout(()=>b.textContent="Share link",1800)}catch{window.prompt("Copy this link:",location.href)}};
let thumbRun=0;
function populateGallery(){const term=$("gallerySearch").value.trim().toLowerCase();const selected=DATES.filter(x=>(x[0]+" "+x[1]+" "+x[2]).toLowerCase().includes(term));const visible=selected.slice(0,visibleDates);$("galleryGrid").innerHTML="";$("galleryCount").textContent=selected.length+" dates";$("showAll").hidden=visible.length>=selected.length;
const rid=++thumbRun;
for(const [date,name,category] of visible){const button=document.createElement("button");button.type="button";button.className="gallery-card";button.innerHTML='<canvas width="180" height="124"></canvas><div class="card-copy"><small></small><strong></strong><span></span></div>';button.querySelector("small").textContent=date;button.querySelector("strong").textContent=name;button.querySelector("span").textContent=category;button.onclick=()=>{$("date").value=date;document.querySelector("#studio").scrollIntoView({behavior:"smooth"});generate(true)};$("galleryGrid").appendChild(button)}
const tasks=[...$("galleryGrid").children];(async()=>{for(let i=0;i<tasks.length;i++){if(rid!==thumbRun)break;try{const [date]=visible[i];const co=coeffs(date),roots=rootsOf(co),view=bounds(roots);const data=await render({co,roots,view,width:180,height:124,palette:"aurora",maxIter:48});if(rid!==thumbRun)break;if(!data.error)paint(tasks[i].querySelector("canvas"),data,roots,view)}catch{}}})()}
$("gallerySearch").oninput=()=>{visibleDates=12;populateGallery()};$("showAll").onclick=()=>{visibleDates=38;populateGallery()};
(async()=>{try{const co=coeffs("25/12/1642"),roots=rootsOf(co),view=bounds(roots);const hero=await render({co,roots,view,width:500,height:500,palette:"aurora",maxIter:80});if(!hero.error)paint($("heroCanvas"),hero,roots,view)}catch{}})();
generate(true);populateGallery();
/* Citation controls are independent of the canvas renderer. */
const citationStrings={
 plain:"Cunha Jr., A. (2026). BirthdayFractal: Newton fractals from calendar dates [Computer software]. GitHub. https://github.com/americocunhajr/BirthdayFractal",
 bibtex:"@misc{CunhaJr2026BirthdayFractal,\n  author       = {Cunha Jr., Americo},\n  title        = {BirthdayFractal: Newton Fractals from Calendar Dates},\n  year         = {2026},\n  howpublished = {Computer software},\n  url          = {https://github.com/americocunhajr/BirthdayFractal}\n}"
};
document.querySelectorAll("[data-copy-citation]").forEach((button)=>{
 button.addEventListener("click",async()=>{
  const value=citationStrings[button.dataset.copyCitation];
  const original=button.innerHTML;
  try{
   if(navigator.clipboard&&window.isSecureContext){await navigator.clipboard.writeText(value);}
   else{
    const area=document.createElement("textarea");area.value=value;area.style.position="fixed";area.style.opacity="0";
    document.body.appendChild(area);area.select();if(!document.execCommand("copy"))throw new Error("Copy unavailable");area.remove();
   }
   button.textContent="Copied!";
  }catch{
   window.prompt("Copy this citation:",value);
   button.textContent="Select and copy";
  }
  setTimeout(()=>{button.innerHTML=original},2000);
 });
});

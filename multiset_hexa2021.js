//
//   Cell Emergent-Ishida Model 
//   2021.10.18
//   2022. 3.31      Copy rigth : Takeshi Ishida 
//
var iterate  = 0  ;       // Iteration numeber 繰り返し数
var cellSize = 8 ;        // Cell spacing (pixel)　セルの間隔（ピクセル）
var meshNum  =100 ;       // Number of vertical and horizontal cells　縦横のセル数
var boxNum   =100000;     // Number of polymerized molecules (number of virtual boxes)　重合分子数（ボックス数）
var molNumMax =100 ;      // Number of molecules polymerized per box　1ボックス当たりの重合する分子数
var w   = 0.650 ;         // Morphological parameter　w 形状パラメータｗ 
var mol = 15 ;            // Number of types of molecules　分子種の数
var x,y ;

var cells   = new Array();   // Cell State　セルの状態
var cells_x = new Array();   // x-coordinate of cell セルのｘ座標
var cells_y = new Array();   // y-coordinate of cell　セルのｙ座標
var cells_boxnum = new Array();   // Box number of polymerized molecules present in each cell　セルに存在する重合分子の番号
var nextCells = new Array(); // For storing the next cell state　次のセル状態の保存
 
var adjacent_cells_x = new Array();   // x sequence number of adjacent cell　隣接セルのｘ配列番号
var adjacent_cells_y = new Array();   // y sequence number of adjacent cell　隣接セルのｙ配列番号
var residualRate = new Array();       // Residual rate by molecular species　各分子の残存率

var box = new Array();                // Polymerized molecule = virtual box 重合分子
var box_num = new Array();            // Number of polymerized molecules of each cell　各セルの重合分子の数
var box_num_pre = new Array();        // Number of polymerized molecules of each cell　(previous time step)　各セルの重合分子の数（1つ前のステップ）
var box_residualRate = new Array();   // Residual rate per polymerized molecule 重合分子の残存率
var box_composition = new Array();    // Molecular configuration within the polymerized molecule 重合分子の構成

var p_entropy = new Array();   // Entropy production エントロピー生成量
var m_entropy = new Array();   // Morphological entropy 形態エントロピー量
var max_m_entropy = new Array();   // Maximum morphological entropy = Configuration entropy when molecules are uniformly distributed 形態エントロピーの最大値（分子が一様に分布しているときのエントロピー）
var total_p_entropy = new Array();   // Time series record of entropy production エントロピー生成量の時系列記録
var total_m_entropy = new Array();   // Time series record of morphological entropy 形態エントロピー量の時系列記録

var canvas;
var ctx;

var buttonStart;
var buttonRandom;
var buttonReset;
var timer1;
var running = false;
 
window.onload = function()
{
    canvas = document.getElementById('Cell Emergent Model');
    ctx = canvas.getContext('2d');
	// Initialization 初期設定
    initCells();
	// Button settings ボタン設定
    buttonStart  = document.getElementById('buttonStart');
    buttonRandom = document.getElementById('buttonRandom');
    buttonReset  = document.getElementById('buttonReset');

    buttonStart.addEventListener('click', onStart, false);
    buttonRandom.addEventListener('click', randomCells, false);
    buttonReset.addEventListener('click', initCells, false);
    canvas.addEventListener('click', canvasClick, false);
};
 
// Processing of "Start" button  開始ボタン
function onStart(){
    if(running){
        clearInterval(timer1);
        buttonStart.value = "Start";
        running = false;
    } else {
		// Iterative calculation 繰り返し計算
        nextGeneration();
		
        timer1 = setInterval("nextGeneration()", 20);
        buttonStart.value = "Stop";
        running = true;
    }
}
 
// Processing of "Initialization" button or ”Reset” button 初期化、リセットボタン
function initCells(){
	var ni, nj ;
	
    ctx.fillStyle = 'rgb( 180, 180, 180)';   // Background color 背景色
    ctx.fillRect(0,0, canvas.width, canvas.height);
	
	// Preparation of array variables 配列変数の準備
    for(i=0;i<=meshNum;i++){
        cells[i]   = new Array();
        cells_x[i] = new Array();
        cells_y[i] = new Array();
        cells_boxnum[i] = new Array();
        nextCells[i] = new Array();
        adjacent_cells_x[i] = new Array();
        adjacent_cells_y[i] = new Array();
		p_entropy[i] = new Array() ;
        box_num[i]       = new Array();
        box_num_pre[i]   = new Array();
	
	       for(j=0;j<=meshNum;j++){
           cells[i][j]            = new Array();
           cells_boxnum[i][j]     = new Array();
           nextCells[i][j]        = new Array();
           adjacent_cells_x[i][j] = new Array();
           adjacent_cells_y[i][j] = new Array();
		   p_entropy[i][j] = new Array() ;
           box_num[i][j]       = new Array();
           box_num_pre[i][j]   = new Array();

	       for(k=0;k<2;k++){
              box_num[i][j][k]       = new Array();
              box_num_pre[i][j][k]   = new Array();
		   }
		}
	}
	
    for(i=1;i<=2;i++){
        box[i]     = new Array();
		box_residualRate[i] = new Array(); 
        box_composition[i]  = new Array();    //
        for(j=0;j<=boxNum;j++){
           box[i][j]     = new Array();
  		   box_residualRate[i][j] = new Array(); 
           box_composition[i][j]  = new Array();    //
	    }
	}
    // Set cell coordinates and initialize the number of cell numerators  セルの座標の設定、セルの分子数の初期化 
    for(i=1;i<=meshNum;i++){
        for(j=1;j<=meshNum;j++){
		 	cells_x[i][j]=(0+((j-1)%2)*cellSize/2+(i-1)*cellSize);
		 	cells_y[i][j]=(0+ (j-1)*cellSize*Math.sqrt(3)/2);
		    p_entropy[i][j] = 0 ;

            for(k=0;k<=mol;k++){
                cells[i][j][k]  = 0;
 			}
        }
    }
	// Initial placement of empty virtual box  空ボックスの初期配置
 	for (i = 1; i <= 2; i++) {
 	for (j = 0; j <= boxNum-1; j++) {
  		box[i][j][0]=0;              // Box state number, 0: empty, 1: polymer state  ボックスの状態番号、０：空、１：高分子状態 
  		box[i][j][1]=(Math.floor(j / meshNum)% meshNum)+1 ;   // x-coordinate of the box ボックスのx座標 
  		box[i][j][2]=(j % meshNum)+1 ;                        // y-coordinate of the box ボックスのy座標 

	}
	}
 	for (i = 1; i <= 2; i++) {
 	    for (j = 0; j <= boxNum-1; j++) {
		 	for (k = 0; k <= 7; k++) {
		 	   box_residualRate[i][j][k] = 0.1;         // Residual rate by orientation in each box  各ボックスのボックスごと、方向別の残存率 
		 	}
	    }
	}  
 	for (i = 1; i <= 2; i++) {
 	    for (j = 0; j <= boxNum-1; j++) {
		 	for (k = 0; k <= molNumMax-1; k++) {
               box_composition[i][j][k]=0;              // Initialization of the composition of each box 各ボックスの組成の初期化
		 	}
		}
	}
	// Relation settings of adjacent cells  隣接セルの設定
    for(i=1;i<=meshNum;i++){
        for(j=1;j<=meshNum;j++){
		 		if ((j-1)%2==1) {
		 		  // 1 :Central 中央
		 		  ni= i;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][0]=ni ;
		 		  adjacent_cells_y[i][j][0]=nj ;

		 		  // 2 :Upper 上
		 		  ni= i;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][1]=ni ;
		 		  adjacent_cells_y[i][j][1]=nj ;

		 		  // 3 :Right 右
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][3]=ni ;
		 		  adjacent_cells_y[i][j][3]=nj ;

		 		  // 4 :Lower 下
		 		  ni= i;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][5]=ni ;
		 		  adjacent_cells_y[i][j][5]=nj ;

		 		  // 5 :Left 左
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][6]=ni ;
		 		  adjacent_cells_y[i][j][6]=nj ;

		 		  // 6 : Upper right 右上
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][2]=ni ;
		 		  adjacent_cells_y[i][j][2]=nj ;

		 		  // 7 : Lower right 右下
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][4]=ni ;
		 		  adjacent_cells_y[i][j][4]=nj ;

		 		  // 8 : Lower left 左下
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;

		 		  // 9 : Upper left 左上
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		}

		 		if ((j-1)%2==0) {
		 		  // 1 : Central 中央
		 		  ni= i;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][0]=ni ;
		 		  adjacent_cells_y[i][j][0]=nj ;

		 		  // 2 : Upper 上
		 		  ni= i;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][2]=ni ;
		 		  adjacent_cells_y[i][j][2]=nj ;

		 		  // 3 : Right 右
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][3]=ni ;
		 		  adjacent_cells_y[i][j][3]=nj ;

		 		  // 4 : Lower 下
		 		  ni= i;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][4]=ni ;
		 		  adjacent_cells_y[i][j][4]=nj ;

		 		  // 5 : Left 左
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][6]=ni ;
		 		  adjacent_cells_y[i][j][6]=nj ;

		 		  // 6 : Upper right 右上
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;

		 		  // 7 : Lower right 右下
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;

		 		  // 8 : Lower left 左下
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][5]=ni ;
		 		  adjacent_cells_y[i][j][5]=nj ;

		 		  // 9 : Upper right 左上
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][1]=ni ;
		 		  adjacent_cells_y[i][j][1]=nj ;
		 		}
        }
    }
	
	// Setting of residual rate of each molecule 各分子の残存率の設定
	residualRate[0] = 0.0;
	residualRate[1] = 0.0;
	residualRate[2] = 0.75; // 0.75
	residualRate[3] = 0.05; // 0.05
	residualRate[4] = residualRate[2];
	residualRate[5] = residualRate[3];
	residualRate[6] = 0.0;
	residualRate[7] = 0.0;
	residualRate[8] = 1.0;
	residualRate[9] = 1.0;
	residualRate[10] = 1.0;
	residualRate[11] = 1.0;
	residualRate[12] = 0.0;
	residualRate[13] = 1.0;
	residualRate[14] = 1.0;
	residualRate[15] = 1.0;

    // Graphic Drawing 描画
    redraw();
}
 
// Processing of "Initialization" button  イニシャルボタン；特定のセル設定
function randomCells(){

    //  Initial molecular configuration 分子の初期配置
    for(x=1;x<=meshNum;x++){
        for(y=1;y<=meshNum;y++){
           cells[x][y][1] = 1000000 ;
           cells[x][y][6] = 100000*0.5;
           cells[x][y][7] = 100000*0.5;
           cells[x][y][12] = 100000 ;
           cells[x][y][13] = 0;
           cells[x][y][14] = 3;  //3 
           cells[x][y][15] = 9;  //9
           cells_boxnum[x][y] = 0;
		   
		   box_num[x][y][1] = 0 ;
		   box_num[x][y][2] = 0 ;
 		   box_num_pre[x][y][1] = 0 ;
		   box_num_pre[x][y][2] = 0 ;
        }
    }
	
	// Initial configuration of molecule 2 and molecule 3　分子２および分子３の初期配置 Case.1
    for(var m=2;m<=3;m++){
    x=50; y=50 ; cells[x][y][m] = 100 ; 
    x=51; y=50 ; cells[x][y][m] = 100 ; 
    x=52; y=50 ; cells[x][y][m] = 100 ; 
    x=53; y=50 ; cells[x][y][m] = 100 ; 
	
    x=50; y=51 ; cells[x][y][m] = 100 ; 
    x=51; y=51 ; cells[x][y][m] = 100 ; 
    x=52; y=51 ; cells[x][y][m] = 100 ; 
    x=53; y=51 ; cells[x][y][m] = 100 ; 
	
    x=49; y=52 ; cells[x][y][m] = 100 ; 
    x=50; y=52 ; cells[x][y][m] = 100 ; 
    x=51; y=52 ; cells[x][y][m] = 100 ; 
    x=52; y=52 ; cells[x][y][m] = 100 ; 
	}

	// Initial configuration of molecule 2 and molecule 3　初期配置 Case.2
/*    for(var m=2;m<=3;m++){
    x=50; y=50 ; cells[x][y][m] = 100 ; 
    x=51; y=50 ; cells[x][y][m] = 100 ; 
	
    x=50; y=51 ; cells[x][y][m] = 100 ; 
    x=51; y=51 ; cells[x][y][m] = 100 ; 
	}
*/
	//  Initial configuration of molecule 2 and molecule 3　初期配置　Case.3
/*    for(var m=2;m<=3;m++){
    x=50; y=50 ; cells[x][y][m] = 100 ; 
    x=51; y=50 ; cells[x][y][m] = 100 ; 
    x=52; y=50 ; cells[x][y][m] = 100 ; 
    x=53; y=50 ; cells[x][y][m] = 100 ; 
	
    x=50; y=51 ; cells[x][y][m] = 100 ; 
    x=51; y=51 ; cells[x][y][m] = 100 ; 
    x=52; y=51 ; cells[x][y][m] = 100 ; 
    x=53; y=51 ; cells[x][y][m] = 100 ; 
	
    x=49; y=52 ; cells[x][y][m] = 100 ; 
    x=50; y=52 ; cells[x][y][m] = 100 ; 
    x=51; y=52 ; cells[x][y][m] = 100 ; 
    x=52; y=52 ; cells[x][y][m] = 100 ; 

    x=49; y=53 ; cells[x][y][m] = 100 ; 
    x=50; y=53 ; cells[x][y][m] = 100 ; 
    x=51; y=53 ; cells[x][y][m] = 100 ; 
    x=52; y=53 ; cells[x][y][m] = 100 ; 
	}
*/	
	// Initial configuration of virtual box ボックスの初期配置
	for (var i = 0; i <= boxNum-1 ; i++) {
	
        if ((box[1][i][1]>45)&&(box[1][i][1]<58)) { 	
        if ((box[1][i][2]>44)&&(box[1][i][2]<58)) { 	
		// Polymerization reaction of informant (polymerized molecule 1)　情報子（重合分子１）の重合反応
	    build_polymer(1,i,6,7,w) ;  
        
		// Polymerization reaction of membrain (polymerized molecule 2)　膜分子（重合分子２）の重合反応
	    //build_polymer(2,i,12,0,1.00) ;  
		}
	    }
	}
	// Counting the number of polymerized molecules in each cell　各セルの重合分子数のカウント
	//calc_polmerNum3(1) ;
	//calc_polmerNum3(2) ;
	
	// Graphic Drawing　描画
    redraw();

}
 
// Redraw process　全体を再描画
function redraw(){
    for(x=1;x<=meshNum;x++){
    for(y=1;y<=meshNum;y++){
	    // Display empty cells　空セルの表示
        //ctx.strokeStyle='rgb(0, 0,0)';
        ctx.strokeStyle='rgb(255,255,255)';
        ctx.lineWidth= 0.5;
	    ctx.fillStyle = 'rgb( 255,255,255)';
        ctx.beginPath();
        ctx.arc( cells_x[x][y]+cellSize,cells_y[x][y]+cellSize , cellSize/2, 0 , 2*Math.PI);
	    ctx.fill();
	    ctx.closePath();
	    ctx.stroke();
	}
	}
	

    for(x=1;x<=meshNum;x++){
    for(y=1;y<=meshNum;y++){
	
	    // Draw of polymerized molecule 1　重合分子1の表示
        if (box_num[x][y][1]>0) drawPolymer(x, y, 1);

	    // Draw of molecule　各分子の表示
  	    // Draw of molecule 1　分子１の表示
        //drawCell(x, y, 1);
	    // Draw of molecule 2　分子２の表示
        //drawCell(x, y, 2);
	    // Draw of molecule 3 分子３の表示
        //drawCell(x, y, 3);
	    // Draw of molecule 4 分子４の表示
        //drawCell(x, y, 4);
	    // Draw of molecule 5　分子５の表示
        //drawCell(x, y, 5);
	    // Draw of molecule 8 分子８の表示
        //drawCell(x, y, 8);
	    // Draw of molecule 9　分子９の表示
        //drawCell(x, y, 9);
	    // Draw of molecule 10　分子10の表示
        //drawCell(x, y, 10);
	    // Draw of molecule 11　分子11の表示
        drawCell(x, y, 11);

	    // Draw of polymerized molecule 2　重合分子2の表示
        if ((box_num[x][y][2]>0)&&(box_num[x][y][2]<=1)) drawPolymer(x, y, 2);
        if ((box_num[x][y][2]>1)&&(box_num[x][y][2]<=4)) drawPolymer(x, y, 3);
        if ((box_num[x][y][2]>4)&&(box_num[x][y][2]<=6)) drawPolymer(x, y, 4);
        if ((box_num[x][y][2]>6))                        drawPolymer(x, y, 5);
		
        // Display of iteration number　繰り返し数の表示
        ctx.fillStyle = 'rgba(225, 225, 225, 0.05)';        // Fill color (translucent)　塗りつぶす色,半透明
        ctx.fillRect(15, meshNum*cellSize*Math.sqrt(3)/2-30, 65, 20);    // Rectangle Drawing　矩形描画
		
		ctx.font = '20pt Arial';
		ctx.fillStyle = 'rgba(0, 0, 0)';
		ctx.fillText(iterate, 20, meshNum*cellSize*Math.sqrt(3)/2-10);
    }
    }
	// Display legend　凡例の表示
/*       if (mm==2) ctx.fillStyle = 'rgb( 255,182,193)';  // Pale pink　薄紅色
       if (mm==3) ctx.fillStyle = 'rgb( 255, 25,147)';    // pink　ピンク色
       if (mm==4) ctx.fillStyle = 'rgb( 255, 25, 25)';    // Yellow 赤色
       if (mm==5) ctx.fillStyle = 'rgb( 225,  0,225)';    // Purple 紫色
*/
	// polymerized molecule 1 重合子１
	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+30 ;
       ctx.fillStyle = 'rgb( 180,225,180)';  // Light green 薄い緑

	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 1", x+15, y+5);	

	// polymerized molecule 2 重合子2
	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+50 ;
       ctx.fillStyle = 'rgb( 255,182,193)';  // Light pink 薄紅色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules in the lattice is less than 1)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+70 ;
       ctx.fillStyle = 'rgb( 255, 25,147)';  // Pink ピンク色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 1 and less than 4)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+90 ;
       ctx.fillStyle = 'rgb( 255, 25, 25)';  // Red 赤色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 4 and less than 6)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+110 ;
       ctx.fillStyle = 'rgb( 225,  0,225)';  // Purple 紫色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 6)", x+15, y+5);	

}
 
// Molecular Drawing 分子の描画
function drawCell(xx, yy, mm){
    if (cells[xx][yy][mm]>0) {
	   //　Color set 色の設定
       if (mm==1) ctx.fillStyle = 'rgb( 255,255,255)'; // White 白
       if (mm==2) ctx.fillStyle = 'rgba(255,224,32,0.5)'; // Yellow 黄色
       if (mm==3) ctx.fillStyle = 'rgba(   0,　32,255,0.5)'; // Blue 青
       if (mm==4) ctx.fillStyle = 'rgba( 255,  0,  0,0.5)'; // Red 赤
       if (mm==5) ctx.fillStyle = 'rgb( 160, 32,255)'; // Purple 紫
       if (mm==6) ctx.fillStyle = 'rgb( 255,208,160)'; // Light pink 薄紅色
       if (mm==7) ctx.fillStyle = 'rgb( 160,128, 60)'; // Brown 茶色
       if (mm==8) ctx.fillStyle = 'rgb(  80,208,255)'; // Light blue 水色
       if (mm==9) ctx.fillStyle = 'rgb(   0,192,  0)'; // Green 緑色
       if (mm==10) ctx.fillStyle = 'rgb(  0,255,  0)'; // Green 緑
       if (mm==11) ctx.fillStyle = 'rgb(255, 96,208)'; // Pink ピンク色
       if (mm==12) ctx.fillStyle = 'rgb(255,224, 32)'; // Yellow 黄色
       
	   ctx.fillStyle = 'rgba(255,224,32,0.8)'; // Yellow 黄色
	   
	   // Cell drawing セルの表示
	   ctx.beginPath();
       //ctx.arc( cells_x[xx][yy]+cellSize, cells_y[xx][yy]+cellSize, cellSize/2*cells[xx][yy][mm]/100, 0 , 2*Math.PI);
       //ctx.arc( cells_x[xx][yy]+cellSize/2*Math.cos(2*3.14/12*mm), cells_y[xx][yy]+cellSize/2*Math.sin(2*3.14/12*mm), cellSize/2*cells[xx][yy][mm]/100, 0 , 2*Math.PI);
       ctx.arc( cells_x[xx][yy]+cellSize,cells_y[xx][yy]+cellSize , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();
    }	
}

// Drawing of polymerized molecules 重合分子の描画
function drawPolymer(xx, yy, mm){
	   //　　Color set 色の設定
       if (mm==1) ctx.fillStyle = 'rgb( 180,225,180)';  // Light green 薄い緑
       if (mm==2) ctx.fillStyle = 'rgb( 255,182,193)';  // Light pink 薄紅色
       if (mm==3) ctx.fillStyle = 'rgb( 255, 25,147)';  // Pink ピンク色
       if (mm==4) ctx.fillStyle = 'rgb( 255, 25, 25)';  // Red 赤色
       if (mm==5) ctx.fillStyle = 'rgb( 225,  0,225)';  // Purple 紫色

       // Cell drawing セルの表示
	   ctx.beginPath();
       ctx.arc( cells_x[xx][yy]+cellSize,cells_y[xx][yy]+cellSize , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();
   
} 

// Reaction, polymerization and diffusion processes 反応・重合・拡散プロセス
function nextGeneration(){
    var i ;
    var molrate1 ;
	var x,y ;
	iterate = iterate +1 ;
	//console.log("iterate = ", iterate) ;

    // Initialization of entropy generation　エントロピー生成量の初期化
    total_p_entropy[iterate] = 0;

	// Molecular Diffusion　分子の拡散
	diff(1) ;
	diff(2) ;
	diff(3) ;
	
	// Molecular Reactions　分子の反応  chem_react(x1,x2,y1,y2,r);  x1 + x2 → y1 + y2  Reaction rate r 反応確率　ｒ
	chem_react(2,0,4,0,5) ; // 5
	chem_react(3,0,5,0,2) ; // 2

	// Molecular Diffusion　分子の拡散
	diff(4) ;
    diff(5) ;
	diff(6) ;
	diff(7) ;
//	diff(10) ;
//	diff(11) ;
	diff(12) ;

	// Molecular Reactions 分子の反応  chem_react(x1,x2,y1,y2,r);  x1 + x2 → y1 + y2  Reaction rate r 反応確率　ｒ
	chem_react(4,1,4,8,100) ;
	chem_react(5,1,5,9,100) ;
	chem_react(4,1,4,11,100) ;

	chem_react(8,0,10,0,-1) ;  // -1；Reacts in the ratio of molecules 6 in the polymerized molecule 重合分子中の分子６の比率で反応する
	chem_react(9,0,10,0,-1) ;  // -1；Reacts in the ratio of molecules 6 in the polymerized molecule 重合分子中の分子６の比率で反応する

	chem_react(8,0,1,0,100) ;
	chem_react(9,0,1,0,100) ;

	chem_react(10,11,1,1,100) ;   // Compare numerators 10 and 11, and the one with the higher number remains. 分子10，11のどちらは多いほうが残る
	
	// Molecular Polymerization 分子の重合
    //var nr   = 6 ;  //8
    //var nr2  = 16 ;  //16

 	for (i = 0; i <= boxNum-1; i++) {
		
		// Polymerization of informant (polymerized molecule 1) 情報子の重合
		if (box[1][i][0]==0) {
		if ((cells[box[1][i][1]][box[1][i][2]][11]>0)) {  // If molecule 11 exists 分子11が存在するなら
		molrate1 = calc_molrate(1,i,6,7) ;   // Calculation of the ratio of molecules in a polymerized molecule 分子の比率の計算
        //molrate1 =w ;
		if (molrate1>0) {
            // Molecular Polymerization 分子の重合
		    //console.log(molrate1);
		    build_polymer(1,i,6,7,molrate1) ;
		}  
		}
		}
	}

	chem_react(14,11,14,1,100) ;
	chem_react(15,11,15,13,100) ;
	chem_react(13,11,11,11,100) ;
	
 	for (i = 0; i <= boxNum-1; i++) {
		//if ((cells[box[2][i][1]][box[2][i][2]][11]<=nr2)&&(cells[box[2][i][1]][box[2][i][2]][11]>=nr)) {  //分子11が一定の範囲で存在するなら
		if ((cells[box[2][i][1]][box[2][i][2]][13]>0)) {  //If molecule 13 exists 分子13が存在するなら
			// Formation of molecules 2 and 3 分子２，３の生成  
	        chem_react2(box[2][i][1],box[2][i][2],11,1,11,2,100) ;
	        chem_react2(box[2][i][1],box[2][i][2],11,1,11,3,100) ;

	        chem_react2(box[2][i][1],box[2][i][2],14,1,14,2,100) ;
	        chem_react2(box[2][i][1],box[2][i][2],14,1,14,3,100) ;

	        chem_react2(box[2][i][1],box[2][i][2],13,1,13,2,100) ;
	        chem_react2(box[2][i][1],box[2][i][2],13,1,13,3,100) ;

       // Polymerization of membrane molecules (polymerized molecule 2)　膜分子（重合分子２）の重合
		if (box[2][i][0]==0) { 
/*	       if (cells[box[2][i][1]][box[2][i][2]][11]>0){
		    console.log("11=",cells[box[2][i][1]][box[2][i][2]][11]);
 	        console.log("13=",cells[box[2][i][1]][box[2][i][2]][13]);
		   }
*/          cells[box[2][i][1]][box[2][i][2]][13]=cells[box[2][i][1]][box[2][i][2]][13]-1;
 		    build_polymer(2,i,12,0,1.00) ;  // Molecular Polymerization　分子の重合
		}
		//}
		//if (box[2][i][0]==1) {  //If molecule 13 exists　分子13が存在するなら

		}
	}
	
	// Counting the number of polymerized molecules in each cell　各セルの重合分子数のカウント
	calc_polmerNum3(1) ;
	calc_polmerNum3(2) ;

	// Diffusion of polymerized molecules　重合分子の拡散
	diff_polymer(1,0.10) ;//0.1
	diff_polymer(2,0.75) ;//0.75
	
	
    // Calculation of entropy in the diffusion process of polymerized molecules 　重合分子の拡散によるエントロピー計算 	
	for (i = 1; i <= meshNum; i++) {
	for (j = 1; j <= meshNum; j++) {
    	for (k = 1; k <= 2; k++) {
    		if (box_num[i][j][k]>=1) { 
			  total_p_entropy[iterate] = total_p_entropy[iterate]   + ((box_num[i][j][k]-box_num_pre[i][j][k])*(box_num[i][j][k]-box_num_pre[i][j][k]))/(box_num[i][j][k]);
			  //total_p_entropy[iterate] = total_p_entropy[iterate]   + ((box_num[i][j][k]-box_num_pre[i][j])*(box_num[i][j][k]-box_num_pre[i][j][k]))/(box_num[i][j][k]);
            }
        }
	}
	}	
	
	// Molecular Removal　分子の除去
	chem_react(4,0,1,0,5) ; // 5
	chem_react(5,0,1,0,5) ; // 5
	//chem_react(8,0,1,0,100) ;
	//chem_react(9,0,1,0,100) ;
	chem_react(10,0,1,0,5) ; // 5

	// Degradation of polymerized molecules　重合分子の分解
 	for (i = 0; i <= boxNum-1; i++) {
		// Degradation of informant (polymerized molecule 1)　情報子（重合分子１）の分解
	    //decomp_polymer(1,i,0.00) ;  // 分子の分解

        // Degradation of membrane molecules (polymerized molecules 2)　膜分子（重合分子２）の分解
		//if ((cells[box[2][i][1]][box[2][i][2]][11]<nr)) {  // 10 //If molecule 11 exists within a certain range　分子11が一定の範囲で存在するなら
		if ((cells[box[2][i][1]][box[2][i][2]][13]<=0)) {  // 10 //If molecule 11 exists within a certain range　分子11が一定の範囲で存在するなら
	    decomp_polymer(2,i,0.05);  // Molecular degradation rate 0.05　分子の分解0.05
		}
	}

	chem_react(11,0,1,0,5) ; // 5
	chem_react(13,0,1,0,75) ;//75 50 100
	
	// Calculation and display of the total number of molecules　全分子数の計算、表示
	//console.log(calc_molNum()) ;

    // Calculation of morphological entropy　形態エントロピーの計算
	for (i = 0; i <= mol; i++) {
        m_entropy[i]=0.0;
        max_m_entropy[i]=0.0;
        }
		
    m_entropy[6] =calc_polmerNum2(1) ;
    m_entropy[12]=calc_polmerNum2(2) ;

	for (i = 1; i <= mol; i++) {
        m_entropy[i]=m_entropy[i] + calc_molNum2(i) ;
        }
	
	for (k = 1; k <= mol; k++) {
		if (m_entropy[k]>0) { 
            max_m_entropy[k]=(m_entropy[k]* Math.log(m_entropy[k]) -m_entropy[k])-(meshNum*meshNum*((m_entropy[k]/(meshNum*meshNum)* Math.log(m_entropy[k]/(meshNum*meshNum)) -m_entropy[k]/(meshNum*meshNum))));
            m_entropy[k]= m_entropy[k]* Math.log(m_entropy[k]) -m_entropy[k];
            }
		}

	for (i = 1; i <= meshNum; i++) {
	for (j = 1; j <= meshNum; j++) {
    	for (k = 1; k <= mol; k++) {
    		if (cells[i][j][k]>=1) { 
    		m_entropy[k]=m_entropy[k] - (cells[i][j][k] * Math.log(cells[i][j][k])-cells[i][j][k]);
            }
        }
	}
	}
	for (k = 1; k <= mol; k++) {
        m_entropy[0]=m_entropy[0] + m_entropy[k];
        max_m_entropy[0]= max_m_entropy[0] + max_m_entropy[k];
    }
		
	total_m_entropy[iterate] = (max_m_entropy[0]- m_entropy[0]);
	//console.log(total_p_entropy[iterate],total_m_entropy[iterate]) ;
	
	// Drawing　描画
    redraw();
		
	// Output entropy value　エントロピーデータの出力
	if (iterate==2000) data_download(iterate) ;   // データを出力する繰り返し数を指定する
}
 
// Canvas click 
function canvasClick(e){
    var xx = e.clientX - canvas.offsetLeft;
    var yy = e.clientY - canvas.offsetTop;
    var col = Math.floor(xx / cellSize)+1 ;
    var row = Math.floor(yy / (cellSize*(Math.sqrt(3)/2)))+1 ;
    if (cells[col][row][1]==0)  cells[col][row][1]=1; else cells[col][row][1]=0;
	drawCell(col, row,1);
}

// Chemical reaction　化学反応
function chem_react(p1,p2,q1,q2,r) {
    var pp1, pp2, qq1,qq2 ;
	var rr ;
	var flg =1;
	
	pp1=0;
	pp2=0;
	qq1=0;
	qq2=0; 
	//rr =r ;
	if (p1>0 ) pp1 =1 ;
	if (p2>0 ) pp2 =1 ;
	if (q1>0 ) qq1 =1 ;
	if (q2>0 ) qq2 =1 ;

  	
	for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
            var s = Math.floor( Math.random()*100);

	        var a = cells[x][y][p1];
 	        if (cells[x][y][p1]>0) {
		    

// 	        if (cells[x][y][p2]>0) {
 	        if (p2>0) {
	           if (cells[x][y][p1]<cells[x][y][p2]) {a=cells[x][y][p1]/(qq1+qq2);} else {a=cells[x][y][p2]/(qq1+qq2); }
	        }	
	        //console.log(" a=",a) ; 
 
            if (r==-1) {
			   rr=calc_molrate2(1,x,y,6,7); 
               //rr=w;
               //console.log("rr=", rr);
			   //if ((rr>0)&&(rr<0.4)) console.log("rr=", rr) ;
			   //console.log(cells[x][y][q1]) ;
			   if (rr>0) {
	           if (p1>0) cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1*rr + a*qq2) ; 
	           if (p1>0) cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1*(1-rr) + a*qq2 ); 
	           //if (p2>0) cells[x][y][p2] = cells[x][y][p2]-a*qq1*r - a*qq2 ; 
	           if (q1>0) cells[x][y][q1] = cells[x][y][q1]+Math.round(a*pp1*rr + a*pp2) ; 
	           if (q1>0) cells[x][y][1] = cells[x][y][1]+Math.round(a*pp1*(1-rr) + a*pp2) ; 
	           //if (q2>0) cells[x][y][q2] = cells[x][y][q2]+a*pp1*r + a*pp2; 
			   //console.log(cells[x][y][q1]) ;
			   
			   // Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
			   if (q1==1) flg=-1; else flg = 1;
			   total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;

			   }
			} else if (s<=r) {

	           if (p1>0) cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1 + a*qq2); 
	           if (p2>0) cells[x][y][p2] = cells[x][y][p2]-Math.round(a*qq1 + a*qq2); 
	           if (q1>0) cells[x][y][q1] = cells[x][y][q1]+Math.round(a*pp1 + a*pp2); 
	           if (q2>0) cells[x][y][q2] = cells[x][y][q2]+Math.round(a*pp1 + a*pp2); 

			   // Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
 			   if (q1==1) flg=-1; else flg = 1;
			   total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;
	        }
			}
		}
	}

}
// Chemical reaction 2 化学反応
function chem_react2(xx,yy,p1,p2,q1,q2,r) {
    var pp1, pp2, qq1,qq2 ;
	var rr ;
	var flg =1;
	
	pp1=0;
	pp2=0;
	qq1=0;
	qq2=0; 
	
	if (p1>0 ) pp1 =1 ;
	if (p2>0 ) pp2 =1 ;
	if (q1>0 ) qq1 =1 ;
	if (q2>0 ) qq2 =1 ;


            var s = Math.floor( Math.random()*100);

	        var a = cells[xx][yy][p1];
 	        if (cells[xx][yy][p1]>0) {
		    
 	        if (p2>0) {
	           if (cells[xx][yy][p1]<cells[xx][yy][p2]) {a=cells[xx][yy][p1]/(qq1+qq2);} else {a=cells[xx][yy][p2]/(qq1+qq2); }
	        }	
 
	        if (p1>0) cells[xx][yy][p1] = cells[xx][yy][p1]-Math.round(a*qq1 + a*qq2); 
	        if (p2>0) cells[xx][yy][p2] = cells[xx][yy][p2]-Math.round(a*qq1 + a*qq2); 
	        if (q1>0) cells[xx][yy][q1] = cells[xx][yy][q1]+Math.round(a*pp1 + a*pp2); 
	        if (q2>0) cells[xx][yy][q2] = cells[xx][yy][q2]+Math.round(a*pp1 + a*pp2); 
	        
			// Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
			if (q1==1) flg=-1; else flg = 1;
			total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;

			}
}

// Molecular Diffusion　分子の拡散
function diff(s) {
	
    for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
	        nextCells[x][y][s] = 0;
		}
	}	
	
	for (var i = 1; i <= meshNum; i++) {
	    for (var j = 1; j <= meshNum; j++) {
            t= cells[i][j][s] ;  
            if (t>0) {
			//  Distribution of molecules to adjacent cells　隣接セルへの分配
            for (var p=0 ; p<=6 ; p++) { 
             	cells[i][j][s] =	cells[i][j][s] - Math.floor(t*(1.0-residualRate[s])/7.0);
 			    nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] = nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] + Math.floor(t*(1.0-residualRate[s])/7.0);
                } 
	  		    // Distribution of remaining molecules to own cell　自身のセルへの残存する分の分配
    		    cells[i][j][s] =cells[i][j][s] - Math.floor(t*residualRate[s]);
    		   	nextCells[i][j][s] = nextCells[i][j][s]+ Math.floor(t*residualRate[s]);
	               // If the divisor has a remainder, it is distributed by a random number.　割り算で余りがでた場合、乱数により割ふりをしてしまう。
     	  	       t= cells[i][j][s] ;  

    		   	   while (t>0) {  
                    p = Math.floor( Math.random()*7);
  	  		        if (p<=6) { 
  		    	    cells[i][j][s] = cells[i][j][s] - 1;
   				    nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] =nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] + 1;

      	            } else {
      		    	cells[i][j][s] = cells[i][j][s] - 1;
      		   	    nextCells[i][j][s] = nextCells[i][j][s]+ 1;
                    }
  	  		       t=t-1 ;
              	  }
			}
        }   //  j_loop
	}       //  i_loop
		 
    for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
			
			// Calculation of entropy production associated with molecular diffusion　分子の拡散に伴うエントロピー生成量の計算
			if (cells[x][y][s]>0) {
                total_p_entropy[iterate] = total_p_entropy[iterate]   + (cells[x][y][s]-nextCells[x][y][s])*(cells[x][y][s]-nextCells[i][j][s])/cells[i][j][s];
			}			
			
	        cells[x][y][s] = nextCells[x][y][s];
		}
	}

}

//   Counting of the ratio of molecules in the polymerized molecule (based on cell coordinates)　高分子中の分子の比率のカウント（セル座標ベース）
function calc_molrate2(boxkind,xx,yy,mol1,mol2) {  
	var n1 ,n2 ;
	var rate ;
	var i,j ;
	n1=0 ;
	n2=0 ;
	
/* 	for (i = 0; i <= boxNum-1; i++) {
        if (box[boxkind][i][0]==1) {
		if ((box[boxkind][i][1]==xx)&&(box[boxkind][i][2]==yy)) {
 	  		for (j = 0; j <= molNumMax-1; j++) {
 	  		   if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	  		   if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 		  	}
 	 	}
		}
	}
*/

	// Measures to speed up the calculation of the ratio of polymerized molecule 1 (calculated from one polymerizer without examining all Boxes)　重合分子１の比率の計算を高速化するための措置（全てのBoxを調べずに一つの重合子から計算）
 	i =cells_boxnum[xx][yy] ;
	if (i>0) {
 	for (j = 0; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	    if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 	}
	}
	
	rate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) rate=0;
	//if ((rate>=1.0)&&(rate<1.4)) console.log(n1,n2,rate) ;
	return rate ;
}

//  Count of the ratio of molecules in the polymerized molecule (based on box number)　高分子中の分子の比率のカウント（ボックス番号ベース）
function calc_molrate(boxkind,boxnumber,mol1,mol2) { 
	var n1 ,n2 ;
	var rate ;
	var xx,yy ;
	var i,j ;

	n1=0 ;
	n2=0 ;
	xx= box[boxkind][boxnumber][1] ;
	yy= box[boxkind][boxnumber][2] ;
	
 	//i =boxnumber ;
/*	for (i = 0; i <= boxNum-1; i++) {
        if ((box[boxkind][i][0]==1)) {
		if ((box[boxkind][i][1]==xx)&&(box[boxkind][i][2]==yy)) {
 	  		for (j = 0; j <= molNumMax-1; j++) {
 	  		   if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	  		   if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 		  	}
 	 	}
		}
	}
*/

	// Measures to speed up the calculation of the ratio of polymerized molecule 1 (calculated from one polymerizer without examining all Boxes)　重合分子１の比率の計算を高速化するための措置（全てのBoxを調べずに一つの重合子から計算）
 	i =cells_boxnum[xx][yy] ;
	if (i>0) {
 	for (j = 0; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	    if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 	}
	}	
	rate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) rate=0;
	//if ((rate>=1.00)&&(rate<1.5)) console.log(xx,yy,n1,n2,rate) ;
	return rate 
}

// Molecular Polymerization　分子の重合
function build_polymer (boxkind,boxnumber,mol1,mol2,rate) {
	var i ;
	var p ;
	var n1 ,n2 ;
	var molrate ;
	n1=0 ;
	n2=0 ;	

 	if (box[boxkind][boxnumber][0]==0)  { //If the box is empty　ボックスが空である

 	if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]>molNumMax) {  // If there is a material molecule　材料となる分子がある
 	if ((mol2==0)||((mol2>0)&&(cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]>molNumMax))) {  // If there is a material molecule　材料となる分子がある
 	
	box[boxkind][boxnumber][0]=1 ;
	
  	for (i = 0; i <= molNumMax-1; i++) {
//      p = Math.floor( Math.random()*100)+1; 	
//      p = Math.floor( Math.random()*100)+1; 	

	  	if ((rate*molNumMax)<=i) {
//	  	if ((rate*100)<=(p)) {
	        box_composition[boxkind][boxnumber][i]=mol2;                    // Description of composition in box　各ボックスの組成の記述 		
		    cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]-1 ;
            n2=n2+1 ;
		} else {
	        box_composition[boxkind][boxnumber][i]=mol1;                    // Description of composition in box　各ボックスの組成の記述 		
		    cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]-1 ;
            n1=n1+1 ;
		}
	}
	// Calculation of entropy production associated with molecular polymerization　分子の重合に伴うエントロピー生成量の計算
	total_p_entropy[iterate] = total_p_entropy[iterate] + molNumMax ;
	
	}
	}
	}
	molrate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) molrate=0;
	if ((molrate>0)&&(molrate<0.4)) console.log(boxkind,boxnumber,n1,n2,molrate,rate) ;

}

// Degradation of polymerized molecules　重合分子の分解
function decomp_polymer (boxkind,boxnumber,decomp_rate) {
	
	var molnum ;
	
 	if (box[boxkind][boxnumber][0]==1)  { //If the box is not empty　ボックスが空でない

    var p = Math.floor( Math.random()*100)+1; 	
	
	if (p<=(decomp_rate*100)) {
	    box[boxkind][boxnumber][0]=0 ;
	
  	    for (var i = 0; i <= molNumMax-1; i++) {
			molnum = box_composition[boxkind][boxnumber][i];
			//if (boxkind==1) console.log("composition=",molnum) ;

	        box_composition[boxkind][boxnumber][i]=0;      // Set the composition of the box to 0 (empty)　各ボックスの組成の記述 		
		    cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][molnum] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][molnum]+1 ;
	    }
    	// Calculation of entropy production associated with the decomposition of polymerized molecules　重合分子の分解に伴うエントロピー生成量の計算
		total_p_entropy[iterate] = total_p_entropy[iterate] + molNumMax ;

	}
	}
}

// Diffusion of polymerized molecules　重合分子の拡散
function diff_polymer(mol,rate) {
	var s ;

	// Measures to speed up the calculation of the ratio of polymerized molecules1　重合分子１の比率の計算を高速化するための措置
	if (mol==1) {	
    for(x=1;x<=meshNum;x++){
        for(y=1;y<=meshNum;y++){
           cells_boxnum[x][y] = 0;
        }
    }	
	}
	
    for (var i = 0; i <= boxNum-1; i++) {
       	//if(box[mol][i][0]==1) {
        
		    var p = Math.floor( Math.random()*100)+1; 	
		    var u = Math.floor( Math.random()*7); 	
	  		
			//s = calc_polmerNum(mol,adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u],adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u]) ;
			s = box_num[adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u]][adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u]][mol];
	  		box_residualRate[mol][i][u] = rate + (1-rate)/0.5 *(1/(1+Math.exp(-1*s))-0.5);

	  		if ((box_residualRate[mol][i][u]*100)<=p) {

	  		    box[mol][i][1]= adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u] ;
		 	    box[mol][i][2]= adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u] ;
	  	        
			
			}
	  	//}
	
	// Measures to speed up the calculation of the ratio of polymerized molecules1　重合分子１の比率の計算を高速化するための措置
	if (mol==1) {	
	if(box[mol][i][0]==1) cells_boxnum[box[mol][i][1]][box[mol][i][2]] = i;
	}

    }  // i_loop
}

// Counting of polymers (per cell)　高分子数のカウント（セル単位）
function calc_polmerNum(boxkind,xx,yy) {  
	var n ;
	var number ;
	n=0 ;
	
 	for (var i = 0; i <= boxNum-1; i++) {
	 	if ((box[boxkind][i][1]==xx)&&(box[boxkind][i][2]==yy)) {
           if (box[boxkind][i][0]==1)  n=n+1 ;
		}
	}
	return n ;
}

// Counting of virtual box (per polymerized molecule)　高分子数のカウント（重合分子単位）
function calc_polmerNum2(boxkind) {  
	var n ;
	var number ;
	n=0 ;
	
 	for (var i = 0; i <= boxNum-1; i++) {
           if (box[boxkind][i][0]==1)  n=n+1 ;
	}
	return n ;
}

// Counting of macromolecules (polymerized molecular units) 高分子数のカウント（重合分子セル単位）
function calc_polmerNum3(boxkind) {  
	var n ;
	var number ;
	n=0 ;
	
	for (var i = 1; i <= meshNum; i++) {
	for (var j = 1; j <= meshNum; j++) {	
 	    box_num_pre[i][j][boxkind] = box_num[i][j][boxkind] ;
        box_num[i][j][boxkind] =0 ;
	}
	}
	
	for (var i = 0; i <= boxNum-1; i++) {
           if (box[boxkind][i][0]==1)  {
			 box_num[box[boxkind][i][1]][box[boxkind][i][2]][boxkind] = box_num[box[boxkind][i][1]][box[boxkind][i][2]][boxkind] +1 ;    
			 n=n+1 ;
	       }
	}
	return n ;
}

// Counting of all molecules　全分子のカウント
function calc_molNum() {  
    var num_particle =0;
		for (var i = 1; i <= meshNum; i++) {
		for (var j = 1; j <= meshNum; j++) {
        for (var s = 1; s <= mol  ; s++) {
       		  num_particle = num_particle +cells[i][j][s];    // Counting of molecules　分子数のカウント
        }
		}   //  j_loop
		}   //  i_loop

	  	for (i = 1; i <= 2; i++) {
	  	for (j = 0; j <= boxNum-1; j++) {
	  	 	if (box[i][j][0]==1){
	 	  		for (var k = 0; k <= molNumMax-1; k++) {
	 	  		  if (box_composition[i][j][k]>0) num_particle = num_particle +1 ;
	 		  	}
	 	 	}
	 	}
		}
    return num_particle ;
}
		
// Counting the number of molecules by molecule type　分子種別の分子数のカウント
function calc_molNum2(molkind) {  
    var num_particle =0;
		for (var i = 1; i <= meshNum; i++) {
		for (var j = 1; j <= meshNum; j++) {
       		  num_particle = num_particle +cells[i][j][molkind];    // Counting of molecules　　分子数のカウント
		}   //  j_loop
		}   //  i_loop

    return num_particle ;
}

// CSV output of entropy data　エントロピーデータのCSV出力
function data_download(data_len) {

    var str = "";      // Create empty string　空の文字列を作成

    for(var i = 1; i<=data_len; i++){
        str += i+","+total_p_entropy[i]+","+total_m_entropy[i]+"\n"; // Generate output data　出力データを作成
    }

    var blob =new Blob([str],{type:"text/csv"}); //Set the above string(str) in the array　配列に上記の文字列(str)を設定
    var link =document.createElement('a');
    link.href = URL.createObjectURL(blob); 
    link.download ="tempdate.csv";
    link.click();
	
}
		
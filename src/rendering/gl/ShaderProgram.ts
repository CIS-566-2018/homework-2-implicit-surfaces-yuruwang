import {vec4, mat4, vec2, vec3} from 'gl-matrix';
import Drawable from './Drawable';
import {gl} from '../../globals';

var activeProgram: WebGLProgram = null;

export class Shader {
  shader: WebGLShader;

  constructor(type: number, source: string) {
    this.shader = gl.createShader(type);
    gl.shaderSource(this.shader, source);
    gl.compileShader(this.shader);

    if (!gl.getShaderParameter(this.shader, gl.COMPILE_STATUS)) {
      throw gl.getShaderInfoLog(this.shader);
    }
  }
};

class ShaderProgram {
  prog: WebGLProgram;

  attrPos: number;

  // unifModel: WebGLUniformLocation;
  // unifModelInvTr: WebGLUniformLocation;
  unifViewInv: WebGLUniformLocation;
  // unifColor: WebGLUniformLocation;
  unifFov: WebGLUniformLocation;
  unifEyePos: WebGLUniformLocation;
  unifTime: WebGLUniformLocation;
  unifDimensions: WebGLUniformLocation;
  unifSampler2D: WebGLUniformLocation;


  constructor(shaders: Array<Shader>) {
    this.prog = gl.createProgram();

    for (let shader of shaders) {
      gl.attachShader(this.prog, shader.shader);
    }
    gl.linkProgram(this.prog);
    if (!gl.getProgramParameter(this.prog, gl.LINK_STATUS)) {
      throw gl.getProgramInfoLog(this.prog);
    }

    // Raymarcher only draws a quad in screen space! No other attributes
    this.attrPos = gl.getAttribLocation(this.prog, "vs_Pos");

    // TODO: add other attributes here
    this.unifFov   = gl.getUniformLocation(this.prog, "u_Fov");
    this.unifDimensions = gl.getUniformLocation(this.prog, "u_Dimensions");
    this.unifTime = gl.getUniformLocation(this.prog, "u_Time");
    this.unifEyePos = gl.getUniformLocation(this.prog, "u_EyePos");
    this.unifViewInv = gl.getUniformLocation(this.prog, "u_ViewInv");
    this.unifSampler2D = gl.getUniformLocation(this.prog, "u_BackgroundTexture");
  }

  use() {
    if (activeProgram !== this.prog) {
      gl.useProgram(this.prog);
      activeProgram = this.prog;
    }
  }

  // TODO: add functions to modify uniforms
  setTexture(t: number) {
    this.use();

    if(this.unifTime != -1)
    {
      gl.uniform1f(this.unifTime, t);
    }
  }

  setTime(t: number) {
    this.use();

    if(this.unifTime != -1)
    {
      gl.uniform1f(this.unifTime, t);
    }
  }

  setDimension(dimensions: vec2) {
    this.use();

    if(this.unifDimensions != -1)
    {
      gl.uniform2fv(this.unifDimensions, dimensions);
    } 
  }

  setFov(fov: number) {
    this.use();

    if(this.unifFov != -1)
    {
      gl.uniform1f(this.unifFov, fov);
    }
  }

  setEyePos(eyePos: vec3) {
    this.use();

    if(this.unifEyePos != -1)
    {
      gl.uniform3fv(this.unifEyePos, eyePos);
    }
  }

  setViewInv(viewMatrix: mat4) {
    this.use();

    if(this.unifViewInv != -1)
    {
      let viewInv: mat4 = mat4.create();
      mat4.invert(viewInv, viewMatrix);
      gl.uniformMatrix4fv(this.unifViewInv, false, viewInv);
    }
  }

  draw(d: Drawable) {
    this.use();
    // if(this.unifSampler2D != -1)
    // {
    //     gl.uniform1i(this.unifSampler2D, /*GL_TEXTURE*/0);
    // }

    if (this.attrPos != -1 && d.bindPos()) {
      gl.enableVertexAttribArray(this.attrPos);
      gl.vertexAttribPointer(this.attrPos, 4, gl.FLOAT, false, 0, 0);
    }

    d.bindIdx();
    gl.drawElements(d.drawMode(), d.elemCount(), gl.UNSIGNED_INT, 0);

    if (this.attrPos != -1) gl.disableVertexAttribArray(this.attrPos);

  }
};

export default ShaderProgram;

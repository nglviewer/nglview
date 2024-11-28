import * as yup from 'yup';

export interface IMessage {
    type: string;
    data?: any;
    ID?: string;
    render_params?: any;
    methodName?: string;
    target?: string;
    args?: any[];
    kwargs?: any;
    component_index?: number;
    repr_index?: number;
    last_child?: boolean;
    reconstruc_color_scheme?: boolean;
    movie_making?: boolean;
    buffers?: any[];
    content?: any;
}

export const IMessageSchema = yup.object().shape({
    type: yup.string().required(),
    data: yup.mixed(),
    ID: yup.string(),
    render_params: yup.mixed(),
    methodName: yup.string(),
    target: yup.string(),
    args: yup.array(),
    kwargs: yup.object(),
    component_index: yup.number(),
    repr_index: yup.number(),
    last_child: yup.boolean(),
    reconstruc_color_scheme: yup.boolean(),
    movie_making: yup.boolean(),
    buffers: yup.array(),
    content: yup.mixed(),
});